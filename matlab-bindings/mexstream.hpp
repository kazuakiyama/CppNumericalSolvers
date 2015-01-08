/* class for redirecting cout to matlab's output stream */
class mexstream : public std::streambuf {
public:
  /* call at beginning of mexFunction() */
  static void install() {
#ifdef __APPLE__
    oldoutbuf(std::cout.rdbuf(&getmexout()));
#endif
  }
  static void restore() {
#ifdef __APPLE__
    std::cout.rdbuf(oldoutbuf());
#endif
  }
 private:
  static mexstream & getmexout() {
    static mexstream mexout;
    return mexout;
  }
  static std::streambuf * oldoutbuf(std::streambuf * oldbuf = NULL) {
    static std::streambuf * oldoutbuf;
    std::streambuf * tmp;
    tmp = oldoutbuf;
    oldoutbuf = oldbuf;
    return tmp;
  }
protected:
  virtual std::streamsize xsputn(const char *s, std::streamsize n) {
    mexPrintf("%.*s",n,s);
    return n;
  }
  virtual int overflow(int c = EOF) {
    if (c != EOF) {
      mexPrintf("%.1s",&c);
    }
    return 1;
  }
};
