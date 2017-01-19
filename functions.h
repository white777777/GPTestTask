#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <limits>
#include <array>
#include <eigen3/Eigen/Dense>

/// Function for regression model
/// \f$ f(t) = e^{-Dt} \f$
/// This is reference function for implementing new functions
class FunctionRef
{
public:
  const static size_t nParams = 1;
  typedef Eigen::Matrix<double, nParams ,1> VParams;
  
private:
  /// \f$ f(t) \f$
  inline double CalcFT (const VParams & params, const double t) const
  {
    const double& D = params[0];
    return exp(-D*t);
  }
  
  /// \f$ \frac{df(t)}{d w_i} \f$ where \f$ w_i \f$ - model param
  double CalcDFDIParam(size_t iParam, const VParams & params, const double t) const
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
  
public:
  /// \f$ \int f(t) dt \f$
  /// This function used only for generating test data and shouldn't be
  /// implemented in another functions
  double calcIntF(const VParams & params, const double t) const
  {
    const double& D = params[0];
    return -exp(-D*t)/D;
  }
  
public:
  // Special for RegressionModelLn
  double calcLnFT (const VParams & params, const double t) const
  {
    const double& D = params[0];
    return -D*t;
  }
  
  // Special for RegressionModelLn
  double calcDFDIParamDivFT(size_t iParam, const VParams & params, const double t) const
  {
    switch(iParam)
    {
      case 0 :
        return -t;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  
  inline double GetParamLowerLimits(size_t iParam) const
  {
    switch(iParam)
    {
      case 0 :
        return 0.0;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }

  inline double GetParamUpperLimits(size_t iParam) const
  {
    return std::numeric_limits<double>::max();
  }
  
  double GetDefaultParam(size_t iParam) const 
  {
    switch(iParam)
    {
      case 0 :
        return 1e-20;
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
  
private:
  /// \f$ f(t) \f$
  inline double CalcFT (const VParams & params, const double t) const
  {
    const double& a = params[0];
    return pow(1+t,-a);
  }
  
  /// \f$ \frac{df(t)}{d w_i} \f$ where \f$ w_i \f$ - model param
  double CalcDFDIParam(size_t iParam, const VParams & params, const double t) const
  {
    switch(iParam)
    {
      case 0 :
        const double& a = params[0];
        return -pow(1+t,-a)*log(-a);
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
public:
  // Special for RegressionModelLn
  double calcLnFT (const VParams & params, const double t) const
  {
    const double& a = params[0];
    return log(pow(1+t, -a));
  }
  
  // Special for RegressionModelLn
  double calcDFDIParamDivFT(size_t iParam, const VParams & params, const double t) const
  {
    switch(iParam)
    {
      case 0 :
        const double& a = params[0];
        return -pow(t+1,-a)*log(t+1);
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  inline double GetParamLowerLimits(size_t iParam) const
  {
    switch(iParam)
    {
      case 0 :
        return 0.0;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  inline double GetParamUpperLimits(size_t iParam) const
  {
    return std::numeric_limits<double>::max();
  }
  
  double GetDefaultParam(size_t iParam) const 
  {
    switch(iParam)
    {
      case 0 :
        return 2;
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
  
private:
  /// \f$ f(t) \f$
  inline double CalcFT (const VParams & params, const double t) const
  {
    const double& D = params[0];
    const double& b = params[1];
    return pow(1+b*D*t, -1.0/b);
  }
  
  /// \f$ \frac{df(t)}{d w_i} \f$ where \f$ w_i \f$ - model param
  double CalcDFDIParam(size_t iParam, const VParams & params, const double t) const
  {
    const double& D = params[0];
    const double& b = params[1];
    switch(iParam)
    {
      case 0 :
        return -t*pow(1+b*D*t, -1.0/b-1);
        break;
      case 1 :
        return pow(1+b*D*t, -1.0/b);
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  public:
    /// \f$ \int f(t) dt \f$
    /// This function used only for generating test data and shouldn't be
    /// implemented in another functions
    double calcIntF(const VParams & params, const double t) const
    {
      const double& D = params[0];
      return -exp(-D*t)/D;
    }
    
  public:
    // Special for RegressionModelLn
    double calcLnFT (const VParams & params, const double t) const
    {
      const double& D = params[0];
      return -D*t;
    }
    
    // Special for RegressionModelLn
    double calcDFDIParamDivFT(size_t iParam, const VParams & params, const double t) const
    {
      switch(iParam)
      {
        case 0 :
          return -t;
          break;
      }
      throw std::invalid_argument("Param index exceeds params count");
    }
    
    
    inline double GetParamLowerLimits(size_t iParam) const
    {
      switch(iParam)
      {
        case 0 :
          return 0.0;
          break;
      }
      throw std::invalid_argument("Param index exceeds params count");
    }
    
    inline double GetParamUpperLimits(size_t iParam) const
    {
      return std::numeric_limits<double>::max();
    }
    
    double GetDefaultParam(size_t iParam) const 
    {
      switch(iParam)
      {
        case 0 :
          return 1e-20;
          break;
      }
      throw std::invalid_argument("Param index exceeds params count");
    }
  };

/// Function for regression model
/// \f$ f(t) = e^{-Dt} \f$
/// This is reference function for implementing new functions
class Function4
{
public:
  const static size_t nParams = 1;
  typedef Eigen::Matrix<double, nParams ,1> VParams;
  
private:
  /// \f$ f(t) \f$
  inline double CalcFT (const VParams & params, const double t) const
  {
    const double& D = params[0];
    return exp(-D*t);
  }
  
  /// \f$ \frac{df(t)}{d w_i} \f$ where \f$ w_i \f$ - model param
  double CalcDFDIParam(size_t iParam, const VParams & params, const double t) const
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
public:
  // Special for RegressionModelLn
  double calcLnFT (const VParams & params, const double t) const
  {
    const double& D = params[0];
    return -D*t;
  }
  
  // Special for RegressionModelLn
  double calcDFDIParamDivFT(size_t iParam, const VParams & params, const double t) const
  {
    switch(iParam)
    {
      case 0 :
        return -t;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  
  inline double GetParamLowerLimits(size_t iParam) const
  {
    switch(iParam)
    {
      case 0 :
        return 0.0;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
  
  inline double GetParamUpperLimits(size_t iParam) const
  {
    return std::numeric_limits<double>::max();
  }
  
  double GetDefaultParam(size_t iParam) const 
  {
    switch(iParam)
    {
      case 0 :
        return 1e-20;
        break;
    }
    throw std::invalid_argument("Param index exceeds params count");
  }
};

//const double& a = params[1];
//return pow((1+t), (-a));
//const double& b = params[2];
//return pow((1+b*D*t), -1/b);
//

#endif // FUNCTIONS_H
