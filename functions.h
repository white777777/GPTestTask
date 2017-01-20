#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <limits>
#include <array>
#include <Eigen/Dense>

/// Function for regression model
/// \f$ f(t) = e^{-Dt} \f$
/// This is reference function for implementing new functions
class FunctionRef
{
public:
  const static size_t nParams = 1;
  typedef Eigen::Matrix<double, nParams ,1> VParams;
  
  /// \f$ f(t) \f$
  inline static double CalcFT (const VParams & params, const double t)
  {
    const double& D = params[0];
    return exp(-D*t);
  }
  
  /// \f$ \frac{df(t)}{d w_i} \f$ where \f$ w_i \f$ - model param
  inline static double CalcDFDIParam(size_t iParam, const VParams & params, const double t)
  {
    switch(iParam)
    {
      case 0 :
        const double& D = params[0];
        return -t*exp(-D*t);
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  /// \f$ \int f(t) dt \f$
  /// This function used only for generating test data and shouldn't be
  /// implemented in another functions
  inline static double calcIntF(const VParams & params, const double t)
  {
    const double& D = params[0];
    return -exp(-D*t)/D;
  }
  
  // Special for RegressionModelLn
  inline static double calcLnFT (const VParams & params, const double t)
  {
    const double& D = params[0];
    return -D*t;
  }
  
  // Special for RegressionModelLn
  inline static double calcDFDIParamDivFT(size_t iParam, const VParams & params, const double t)
  {
    switch(iParam)
    {
      case 0 :
        return -t;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  
  inline static double GetParamLowerLimits(size_t iParam)
  {
    switch(iParam)
    {
      case 0 :
        return 0.0;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }

  inline static double GetParamUpperLimits(size_t iParam)
  {
    return std::numeric_limits<double>::max();
  }
  
  inline static double GetDefaultParam(size_t iParam) 
  {
    switch(iParam)
    {
      case 0 :
        return 0.0001;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
};

typedef FunctionRef Function1;

/// Function for regression model
/// \f$ f(t) = (1+t)^{-a} \f$
/// This is reference function for implementing new functions
class Function2
{
public:
  const static size_t nParams = 1;
  typedef Eigen::Matrix<double, nParams ,1> VParams;
  
  /// \f$ f(t) \f$
  inline static double CalcFT (const VParams & params, const double t)
  {
    const double& a = params[0];
    return pow(1+t,-a);
  }
  
  /// \f$ \frac{df(t)}{d w_i} \f$ where \f$ w_i \f$ - model param
  static double CalcDFDIParam(size_t iParam, const VParams & params, const double t)
  {
    const double& a = params[0];
    switch(iParam)
    {
      case 0 :
        return -pow(1+t,-a)*log(t+1);
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  inline static double GetParamLowerLimits(size_t iParam)
  {
    switch(iParam)
    {
      case 0 :
        return 1e-300;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  inline static double GetParamUpperLimits(size_t iParam)
  {
    return std::numeric_limits<double>::max();
  }
  
  inline static double GetDefaultParam(size_t iParam) 
  {
    switch(iParam)
    {
      case 0 :
        return 0.01;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
};

/// Function for regression model
/// \f$ f(t) = (1+b*D*t)^{-1/b} \f$
/// This is reference function for implementing new functions
class Function3
{
public:
  const static size_t nParams = 2;
  typedef Eigen::Matrix<double, nParams ,1> VParams;
  
  /// \f$ f(t) \f$
  inline static double CalcFT (const VParams & params, const double t)
  {
    const double& D = params[0];
    const double& b = params[1];
    return pow(1+b*D*t, -1.0/b);
  }
  
  /// \f$ \frac{df(t)}{d w_i} \f$ where \f$ w_i \f$ - model param
  inline static double CalcDFDIParam(size_t iParam, const VParams & params, const double t)
  {
    const double& D = params[0];
    const double& b = params[1];
    const double dbt = D*b*t;
    const double dbtp1 = dbt+1;
    switch(iParam)
    {
      case 0 :
        return -t*pow(dbtp1, -1.0/b-1);
        break;
      case 1 :
        return 1/(b*b) *(-dbt+dbtp1*log(dbtp1))*pow(dbtp1, -(b+1)/b);
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }

  inline static double GetParamLowerLimits(size_t iParam)
  {
    switch(iParam)
    {
      case 0 :
        return 0.0;
        break;
      case 1:
        return 0.0;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  inline static double GetParamUpperLimits(size_t iParam)
  {
    return std::numeric_limits<double>::max();
  }
  
  inline static double GetDefaultParam(size_t iParam) 
  {
    switch(iParam)
    {
      case 0 :
        return 1e-5;
        break;
      case 1 :
        return 2.0;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
};

/// Function for regression model
/// \f$ f(t) = (1+t)^{-a}, t<tau \f$
/// \f$ f(t) = (1+b*D(t-tau))^{-1/b}, t>=tau \f$
/// This is reference function for implementing new functions
class Function4
{
public:
  const static size_t nParams = 3;
  typedef Eigen::Matrix<double, nParams ,1> VParams;
  
  /// \f$ f(t) \f$
  inline static double CalcFT (const VParams & params, const double t)
  {
    const double& b = params[0];
    const double& a = params[1];
    const double& tau = params[2];
    const double D=a*pow(t+1,-a-1);

    if(t<tau)
    {
      Function2::VParams tmp;
      tmp[0] = a;
      return Function2::CalcFT(tmp, t);
    }
    else
    {
      return Function3::CalcFT({D, b}, t-tau);
    }
  }
  
  /// \f$ \frac{df(t)}{d w_i} \f$ where \f$ w_i \f$ - model param
  inline static double CalcDFDIParam(size_t iParam, const VParams & params, const double t)
  {
    const double& b = params[0];
    const double& a = params[1];
    const double& tau = params[2];
    const double D=a*pow(t+1,-a-1);
    
    switch(iParam)
    {
      case 0 :
        if(t<tau)
          return 0;
        else
          return Function3::CalcDFDIParam(0, {D, b}, t-tau);
        break;
      case 1:
        if(t<tau)
        {
          Function2::VParams tmp;
          tmp[0] = a;
          return Function2::CalcDFDIParam(0, tmp, t);
        }
        else
          return 0;
        break;
      case 2:
        if(t<tau)
          return 0;
        else
          return D*pow(D*b*(t-tau)+1, -(b+1)/b);
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }

  inline static double GetParamLowerLimits(size_t iParam)
  {
    switch(iParam)
    {
      case 0 :
      case 1 :
      case 2 :
        return 0.0;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  inline static double GetParamUpperLimits(size_t iParam)
  {
    return std::numeric_limits<double>::max();
  }
  
  inline static double GetDefaultParam(size_t iParam) 
  {
    switch(iParam)
    {
      case 0 :
        return 2.0;
        break;
      case 1 :
        return 0.01;
        break;
      case 2 :
        return 10000;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
};


#endif // FUNCTIONS_H
