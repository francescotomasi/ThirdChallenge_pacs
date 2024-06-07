#ifndef FUN_MUPARSER
#define FUN_MUPARSER

#include <muParser.h>

#include <vector>
#include <memory>
#include <string>
#include <cmath>

class MuparserFun
{
public:
  MuparserFun(const MuparserFun &m)
    : m_parser(m.m_parser)
  {
    if (m.m_var.size() !=2)
      {
        std::cerr << "Error: MuparserFun copy constructor only works for 2 variables" << std::endl;
      }
    m_var.resize(m.m_var.size());
    m_parser.DefineVar("x", &m_var[0]);
    m_parser.DefineVar("y", &m_var[1]);
    m_parser.DefineConst("pi", M_PI);

  };

  MuparserFun(const std::string &s, size_t size)
  {
    if (size != 2)
      {
        std::cerr << "Error: MuparserFun constructor only works for 2 variables" << std::endl;
      }
    m_var.resize(size);
    try
      {
        m_parser.DefineVar("x", &m_var[0]);
        m_parser.DefineVar("y", &m_var[1]);
        m_parser.DefineConst("pi", M_PI);
        m_parser.SetExpr(s);
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
  };

  double
  operator()(const double &x, const double &y)
  {
    m_var[0] = x;
    m_var[1] = y;
    double z = m_parser.Eval();
    return z;
  };

private:
  std::vector<double> m_var;
  mu::Parser m_parser;
};

#endif // FUN_MUPARSER