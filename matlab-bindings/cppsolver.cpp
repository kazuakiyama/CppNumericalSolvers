// Patrick Wieschollek
// for compiling download eigen and call the m-file "make.m" inside Matlab
#include <iostream>
#include <vector>
#include <memory>
#include <Eigen/Eigen>
#include "CppNumericalSolvers/src/CppNumericalSolvers.hpp"
#include "mex.h"
#include "mexstream.hpp"

/* usage: [x,fx] = cppsolver(x0,@objective,[args])
args = 'gradient', @gradient
       'solver', ["gradientdescent"|"newton"|"cg"|"asa_cg"|"bfgs"|"l-bfgs"|"l-bfgs-b"]
       'skip_gradient_check', [default:"false"|"true"]
       'skip_hessian_check', [default:"true"|"false"]
*/

using namespace pwie;

solver_id selected_solver;
char *objective_name;
char *gradient_name;
char *hessian_name;
size_t in_rows, in_cols;
bool has_gradient;
bool has_hessian;
bool has_upperbound;
bool has_lowerbound;
Eigen::VectorXd solution, gradient;
Eigen::VectorXd upper, lower;

struct MexFunctor {
  typedef double Scalar;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = 1
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,ValuesAtCompileTime> JacobianType;
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,InputsAtCompileTime> HessianType;

  // create function
  double f(const Eigen::Matrix<double,InputsAtCompileTime,1> & x) const {
    mxArray * objective_ans, *objective_param[1];
    objective_param[0] = mxCreateDoubleMatrix(x.rows(), x.cols(), mxREAL);
    const double *constVariablePtr = &x(0);
    memcpy(mxGetPr(objective_param[0]), constVariablePtr,
           mxGetM(objective_param[0]) * mxGetN(objective_param[0]) * sizeof(*constVariablePtr));
    mexCallMATLAB(1, &objective_ans, 1, objective_param, objective_name) ;
    return mxGetScalar(objective_ans);
  }

  // create derivative of function
  void gradient(const InputType & x, JacobianType & grad) const {
    if (has_gradient) {
      mxArray * objective_ans, *objective_param[1];
      objective_param[0] = mxCreateDoubleMatrix(x.rows(), x.cols(), mxREAL);
      memcpy(mxGetPr(objective_param[0]), x.data(),
             mxGetM(objective_param[0]) * mxGetN(objective_param[0]) * sizeof(double));

      mexCallMATLAB(1, &objective_ans, 1, objective_param, gradient_name) ;
      size_t r = mxGetM(objective_ans);
      size_t c = mxGetN(objective_ans);
      if ((in_rows != r) || (1 != c)) {
        char error_msg[256];
        sprintf(error_msg, "Wrong format of gradient! The correct format is %zu x %zu,"
                " but %zu x %zu was given", in_rows, in_cols, r, c);
        mexErrMsgIdAndTxt("MATLAB:cppsolver", error_msg);
      }

      grad = Eigen::Map<Eigen::VectorXd>(mxGetPr(objective_ans), mxGetM(objective_ans) );
      //mexPrintf("gradient: %f %f\n",gradient[0],gradient[1]);
    } else
      computeGradient<MexFunctor>(*this, x, grad);
  }

  void hessian(const InputType & x, HessianType & hes) const {
    if (has_hessian) {
      // use provided hessian
      mxArray * objective_ans, *objective_param[1];
      objective_param[0] = mxCreateDoubleMatrix(x.rows(), x.cols(), mxREAL);
      memcpy(mxGetPr(objective_param[0]), x.data(),
             mxGetM(objective_param[0]) * mxGetN(objective_param[0]) * sizeof(double));

      mexCallMATLAB(1, &objective_ans, 1, objective_param, hessian_name) ;
      size_t r = mxGetM(objective_ans);
      size_t c = mxGetN(objective_ans);
      if ((in_rows != r) || (in_rows != c)) {
        char error_msg[256];
        sprintf(error_msg, "Wrong format of hessian! The correct format is %zu x %zu, but %zu x %zu was given",
                in_rows, in_rows, r, c);
        mexErrMsgIdAndTxt("MATLAB:cppsolver", error_msg);
      }

      hes = Eigen::Map<Eigen::MatrixXd>(mxGetPr(objective_ans), mxGetM(objective_ans), mxGetN(objective_ans));
      //mexPrintf("hessian: %f %f\n",gradient[0],gradient[1]);
    } else {
      // numerical approximation of hessian
      hes = Eigen::MatrixXd::Zero(in_rows, in_rows);
      computeHessian<MexFunctor>(*this, x, hes);
    }
  }
};

void mexFunction(int outLen, mxArray *outArr[], int inLen, const mxArray *inArr[])
{
  mexstream::install();

  has_gradient             = false;
  has_hessian              = false;
  bool skip_gradient_check = false;
  bool skip_hessian_check  = true;
  has_upperbound           = false;
  has_lowerbound           = false;
  selected_solver          = SOLVER_BFGS;

  if (inLen < 2) {
    mexErrMsgIdAndTxt("MATLAB:cppsolver", "this function need at leat one parameter");
  }

  // PARSING PARAMETERS

  // initial solution
  in_rows = mxGetM(inArr[0]);
  in_cols = mxGetN(inArr[0]);
  if (in_cols < 1 || in_rows == 0) {
    char error_msg[256];
    sprintf(error_msg, "The first argument has to be the inital guess x0 (format: n x 1),"
            " but the input format is %zu x %zu", in_rows, in_cols);
    mexErrMsgIdAndTxt("MATLAB:cppsolver", error_msg);
  }

  solution = Eigen::Map<Eigen::VectorXd>(mxGetPr(inArr[0]), in_rows * in_cols);

  // function handle "@objective"
  if (mxGetClassID(inArr[1]) != mxFUNCTION_CLASS) {
    mexErrMsgIdAndTxt("MATLAB:cppsolver", "the second arguments has to be the handle of the function (@objective)");
  }

  mxArray *objective_ans, *objective_param[1];

  // get name of objective
  objective_param[0] = const_cast<mxArray *>( inArr[1] );
  mexCallMATLAB(1, &objective_ans, 1, objective_param, "char") ;
  objective_name = mxArrayToString(objective_ans);

  //mexPrintf("Found objective function: %s\n", objective_name);

  // there are some parameters
  if ((inLen % 2) != 0) {
    mexErrMsgIdAndTxt("MATLAB:cppsolver", "optional arguments have to be passed by 'key','value'.");
  }

  for (int arg = 2; arg < inLen; arg += 2) {
    if (!mxIsChar(inArr[arg])) {
      mexErrMsgIdAndTxt("MATLAB:cppsolver", "optional argument keys have to be strings");
    }
    char *key_str = mxArrayToString(inArr[arg]);
    //printf("parsing key: %s\n", key_str);

    if (strcmp(key_str, "gradient") == 0) {
      // extract gradient name
      if (mxGetClassID(inArr[arg + 1]) != mxFUNCTION_CLASS) {
        mexErrMsgIdAndTxt("MATLAB:cppsolver", "the argument following 'gradient' has to a function handle (@gradient)");
      }
      objective_param[0] = const_cast<mxArray *>( inArr[arg + 1] );
      mexCallMATLAB(1, &objective_ans, 1, objective_param, "char") ;
      gradient_name =   mxArrayToString(objective_ans);
      has_gradient = true;
      //mexPrintf("Found gradient function: %s\n", gradient_name);
    }
    else if (strcmp(key_str, "solver") == 0) {
      if (!mxIsChar( inArr[arg + 1])) {
        mexErrMsgIdAndTxt("MATLAB:cppsolver", "solver name has to be a string");
      }
      char *solver_str = mxArrayToString(inArr[arg + 1]);
      if (strcmp(solver_str, "gradientdescent") == 0) {
        selected_solver = SOLVER_GRADD;
      } else if (strcmp(solver_str, "cg") == 0) {
        selected_solver = SOLVER_CG;
      } else if (strcmp(solver_str, "asa_cg") == 0) {
        selected_solver = SOLVER_ASA_CG;
      } else if (strcmp(solver_str, "bfgs") == 0) {
        selected_solver = SOLVER_BFGS;
      } else if (strcmp(solver_str, "l-bfgs") == 0) {
        selected_solver = SOLVER_LBFGS;
      } else if (strcmp(solver_str, "l-bfgs-b") == 0) {
        selected_solver = SOLVER_LBFGSB;
      } else if (strcmp(solver_str, "newton") == 0) {
        selected_solver = SOLVER_NEWTON;
      } else {
        char error_msg[256];
        sprintf(error_msg, "unknown solver %s", solver_str);
        mexErrMsgIdAndTxt("MATLAB:cppsolver", error_msg);
      }
    }
    else if (strcmp(key_str, "skip_gradient_check") == 0) {
      if (!mxIsChar( inArr[arg + 1])) {
        mexErrMsgIdAndTxt("MATLAB:cppsolver", "the value of the key 'skip_gradient_check' has to be a string");
      }
      char *txt = mxArrayToString(inArr[arg + 1]);
      if (strcmp(txt, "true") == 0) {
        skip_gradient_check = true;
      }
    }
    else if (strcmp(key_str, "skip_hessian_check") == 0) {
      if (!mxIsChar( inArr[arg + 1])) {
        mexErrMsgIdAndTxt("MATLAB:cppsolver", "the value of the key 'skip_hessian_check' has to be a string");
      }
      char *txt = mxArrayToString(inArr[arg + 1]);
      if (strcmp(txt, "false") == 0) {
        skip_hessian_check = false;
      }
    }
    else if (strcmp(key_str, "hessian") == 0) {
      // extract hessian name
      if (mxGetClassID(inArr[arg + 1]) != mxFUNCTION_CLASS) {
        mexErrMsgIdAndTxt("MATLAB:cppsolver", "the argument following 'hessian' has to a function handle (@hessian)");
      }
      objective_param[0] = const_cast<mxArray *>( inArr[arg + 1] );
      mexCallMATLAB(1, &objective_ans, 1, objective_param, "char") ;
      hessian_name =   mxArrayToString(objective_ans);
      has_hessian = true;
      //mexPrintf("Found hessian function: %s\n", hessian_name);
    }
    else if (strcmp(key_str, "ub") == 0 || strcmp(key_str, "lb") == 0) {
      // extract UpperBound
      size_t b_in_rows = mxGetM(inArr[arg + 1]);
      size_t b_in_cols = mxGetN(inArr[arg + 1]);

      if ((b_in_cols != 1) || (b_in_rows != in_rows)) {
        char error_msg[256];
        sprintf(error_msg, "The format of the bounds '%s' have to match the format of "
                "the inital guess x0 (format: %zu x 1), but the input format is %zu x %zu",
                key_str, in_rows, b_in_rows, b_in_cols);
        mexErrMsgIdAndTxt("MATLAB:cppsolver", error_msg);
      }

      auto tmp = Eigen::Map<Eigen::VectorXd>(mxGetPr(inArr[arg + 1]),
                                             mxGetM(inArr[arg + 1]) * mxGetN(inArr[arg + 1]));

      if (strcmp(key_str, "lb") == 0) {
        lower = tmp;
        has_lowerbound = true;
      }
      else {
        upper = tmp;
        has_upperbound = true;
      }
    }
    else {
      char error_msg[256];
      sprintf(error_msg, "unknown argument %s", key_str);
      mexErrMsgIdAndTxt("MATLAB:cppsolver", error_msg);
    }
  }

  std::unique_ptr<ISolver<MexFunctor> > g(getSolver<MexFunctor>(selected_solver));

  // check gradient
  if (has_gradient && !skip_gradient_check) {
    Eigen::VectorXd dx = solution;
    g->gradient(solution, dx);
    if (!checkGradient<MexFunctor>(*g, solution, dx)) {
      mexErrMsgIdAndTxt("MATLAB:cppsolver:gradient_check",
                        "your gradient seems to be not correct! You can skip this test by using "
                        "the arguments \"'skip_gradient_check','true'\"");
    }
  }
  // check hessian
  if (!skip_hessian_check) {
    Eigen::MatrixXd hes = Eigen::MatrixXd::Zero(in_rows, in_rows);
    Eigen::MatrixXd hes2 = Eigen::MatrixXd::Zero(in_rows, in_rows);
    g->hessian(solution, hes);
    computeHessian<MexFunctor>(*g, solution, hes2);
    const double diff = static_cast<Eigen::MatrixXd>(hes - hes2).norm() ;
    if (diff > 1e-3) {
      char error_msg[256];
      sprintf(error_msg, "Your hessian is probably not correct or the objective function is "
              "obnoxious(diff to numerical approx: %f)! You can skip this test by removing the "
              "arguments \"'skip_hessian_check','false'\"", diff);
      mexErrMsgIdAndTxt("MATLAB:cppsolver", error_msg);
    }
  }

  gradient = solution;

  if (has_lowerbound)
    g->setLowerBound(lower);

  if (has_upperbound)
    g->setUpperBound(upper);

  g->solve(solution);

  // prepare solution
  outArr[0] = mxCreateDoubleMatrix(solution.rows(), solution.cols(), mxREAL);
  double *constVariablePtr = &solution(0);
  memcpy(mxGetPr(outArr[0]), constVariablePtr, mxGetM(outArr[0]) * mxGetN(outArr[0]) * sizeof(*constVariablePtr));
  outArr[1] = mxCreateDoubleScalar(g->f(solution));

}

