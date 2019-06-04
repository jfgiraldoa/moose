#pragma once

#include "SinglePhaseFluidProperties.h"

class FlinakFluidProperties;

template <>
InputParameters validParams<FlinakFluidProperties>();

/**
 * Fluid properties for 0.465 LiF - 0.115 NaF - 0.42 KF (flinak) \cite richard.
 */
class FlinakFluidProperties : public SinglePhaseFluidProperties
{
public:
  FlinakFluidProperties(const InputParameters & parameters);

  /**
   * Fluid name
   *
   * @return "flinak"
   */
  virtual std::string fluidName() const override;

  /**
   * Pressure from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @return pressure (Pa)
   */
  virtual Real p_from_v_e(Real v, Real e) const override;

  /**
   * Pressure and its derivatives from specific volume and specific internal energy
   *
   * @param[in] v        specific volume (m$^3$/kg)
   * @param[in] e        specific internal energy (J/kg)
   * @param[out] p       pressure (Pa)
   * @param[out] dp_dv   derivative of pressure w.r.t. specific volume
   * @param[out] dp_de   derivative of pressure w.r.t. specific internal energy
   */
  virtual void p_from_v_e(Real v, Real e, Real & p, Real & dp_dv, Real & dp_de) const override;

  /**
   * Temperature from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @return temperature (K)
   */
  virtual Real T_from_v_e(Real v, Real e) const override;

  /**
   * Temperature and its derivatives from specific volume and specific internal energy
   *
   * @param[in] v        specific volume (m$^3$/kg)
   * @param[in] e        specific internal energy (J/kg)
   * @param[out] T       temperature (K)
   * @param[out] dT_dv   derivative of temperature w.r.t. specific volume
   * @param[out] dT_de   derivative of temperature w.r.t. specific internal energy
   */
  virtual void T_from_v_e(Real v, Real e, Real & T, Real & dT_dv, Real & dT_de) const override;

  using SinglePhaseFluidProperties::cp_from_v_e;

  /**
   * Isobaric specific heat from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @return isobaric specific heat (J/kg.K)
   */
  virtual Real cp_from_v_e(Real v, Real e) const override;

  using SinglePhaseFluidProperties::cv_from_v_e;

  /**
   * Isochoric specific heat from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @return isochoric specific heat (J/kg.K)
   */
  virtual Real cv_from_v_e(Real v, Real e) const override;

  using SinglePhaseFluidProperties::mu_from_v_e;

  /**
   * Dynamic viscosity from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @return dynamic viscosity (Pa.s)
   */
  virtual Real mu_from_v_e(Real v, Real e) const override;

  using SinglePhaseFluidProperties::k_from_v_e;

  /**
   * Thermal conductivity from specific volume and specific internal energy
   *
   * @param[in] v   specific volume (m$^3$/kg)
   * @param[in] e   specific internal energy (J/kg)
   * @return thermal conductivity (W/m.K)
   */
  virtual Real k_from_v_e(Real v, Real e) const override;

  /**
   * Density from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @return density (kg/m$^3$)
   */
  virtual Real rho_from_p_T(Real p, Real T) const override;

  /**
   * Density and its derivatives from pressure and temperature
   *
   * @param[in] p          pressure (Pa)
   * @param[in] T          temperature (K)
   * @param[out] rho       density (kg/m$^3$)
   * @param[out] drho_dp   derivative of density w.r.t. pressure
   * @param[out] drho_dT   derivative of density w.r.t. temperature
   */
  virtual void
  rho_from_p_T(Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const override;

  /**
   * Specific volume from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @return specific volume (m$^3$/kg)
   */
  virtual Real v_from_p_T(Real p, Real T) const override;

  virtual DualReal v_from_p_T(const DualReal & p, const DualReal & T) const override;

  /**
   * Specific volume and its derivatives from pressure and temperature
   *
   * @param[in] p          pressure (Pa)
   * @param[in] T          temperature (K)
   * @param[out] v         specific volume (m$^3$/kg)
   * @param[out] dv_dp     derivative of specific volume w.r.t. pressure
   * @param[out] dv_dT     derivative of specific volume w.r.t. temperature
   */
  virtual void v_from_p_T(Real p, Real T, Real & v, Real & dv_dp, Real & dv_dT) const override;

  /**
   * Specific enthalpy from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @return specific enthalpy (J/kg)
   */
  virtual Real h_from_p_T(Real p, Real T) const override;

  /**
   * Specific enthalpy and its derivatives from pressure and temperature
   *
   * @param[in] p        pressure (Pa)
   * @param[in] T        temperature (K)
   * @param[out] h       specific enthalpy (J/kg)
   * @param[out] dh_dp   derivative of specific enthalpy w.r.t. pressure
   * @param[out] dh_dT   derivative of specific enthalpy w.r.t. temperature
   */
  virtual void h_from_p_T(Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const override;

  /**
   * Specific internal energy from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @return specific internal energy (J/kg)
   */
  virtual Real e_from_p_T(Real p, Real T) const override;

  /**
   * Specific internal energy and its derivatives from pressure and temperature
   *
   * @param[in] p        pressure (Pa)
   * @param[in] T        temperature (K)
   * @param[out] e       specific internal energy (J/kg)
   * @param[out] de_dp   derivative of specific internal energy w.r.t. pressure
   * @param[out] de_dT   derivative of specific internal energy w.r.t. temperature
   */
  virtual void e_from_p_T(Real p, Real T, Real & e, Real & de_dp, Real & de_dT) const override;

  using SinglePhaseFluidProperties::beta_from_p_T;

  /**
   * Thermal expansion coefficient from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @return thermal expansion coefficient (1/K)
   */
  virtual Real beta_from_p_T(Real p, Real T) const override;

  /**
   * Isobaric specific heat capacity from pressure and temperature
   *
   * @param p   pressure (Pa)
   * @param T   temperature (K)
   * @return isobaric specific heat (J/kg/.K)
   */
  virtual Real cp_from_p_T(Real p, Real T) const override;

  /**
   * Isobaric specific heat capacity and its derivatives from pressure and temperature
   *
   * @param[in] p       pressure (Pa)
   * @param[in] T       temperature (K)
   * @param[out] cp     isobaric specific heat (J/kg/K)
   * @param[out] dcp_dp derivative of isobaric specific heat w.r.t. pressure (J/kg/K/Pa)
   */
  virtual void cp_from_p_T(Real p, Real T, Real & cp, Real & dcp_dp, Real & dcp_dT) const override;

  using SinglePhaseFluidProperties::cv_from_p_T;

  /**
   * Isochoric specific heat capacity from pressure and temperature
   *
   * @param p   pressure (Pa)
   * @param T   temperature (K)
   * @return isochoric specific heat (J/kg.K)
   */
  virtual Real cv_from_p_T(Real p, Real T) const override;
  virtual void cv_from_p_T(Real p, Real T, Real & cv, Real & dcv_dp, Real & dcv_dT) const override;

  /**
   * Molar mass
   *
   * @return molar mass (kg/mol)
   */
  virtual Real molarMass() const override;

  /**
   * Thermal conductivity from pressure and temperature
   *
   * @param p   pressure (Pa)
   * @param T   temperature (K)
   * @return thermal conductivity  (W/m.K)
   */
  virtual Real k_from_p_T(Real p, Real T) const override;

  /**
   * Thermal conductivity and its derivatives wrt pressure and temperature
   *
   * @param p           pressure (Pa)
   * @param T           temperature (K)
   * @param[out]  k     thermal conductivity  (W/m.K)
   * @param[out]  dk_dp derivative of thermal conductivity wrt pressure
   * @param[out]  dk_dT derivative of thermal conductivity wrt temperature
   */
  virtual void k_from_p_T(Real p, Real T, Real & k, Real & dk_dp, Real & dk_dT) const override;

  /**
   * Dynamic viscosity from pressure and temperature
   *
   * @param p   pressure (Pa)
   * @param T   temperature (K)
   * @return dynamic viscosity (Pa.s)
   */
  virtual Real mu_from_p_T(Real p, Real T) const override;

  /**
   * Dynamic viscosity and its derivatives wrt pressure and temperature
   *
   * @param p             pressure (Pa)
   * @param T             temperature (K)
   * @param[out] mu       viscosity (Pa.s)
   * @param[out] dmu_dp   derivative of viscosity wrt pressure
   * @param[out] dmu_dT   derivative of viscosity wrt temperature
   */
  virtual void
  mu_from_p_T(Real p, Real T, Real & mu, Real & dmu_drho, Real & dmu_dT) const override;

protected:
  /// Derivative of density with respect to pressure at fixed temperature
  const Real & _drho_dp;

  /// Derivative of density with respect to temperature at fixed pressure
  const Real _drho_dT;

  /// Atmospheric pressure, Pa
  const Real _p_atm;

  /// specific heat at constant pressure
  const Real _cp;

  /// additive constant to rho(P, T) correlation
  const Real _c0;

  /// derivative of pressure with respect to temperature at constant specific volume
  const Real _dp_dT_at_constant_v;
};
