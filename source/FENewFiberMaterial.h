#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! A material class describing the active fiber contraction
class FENewActiveFiberContraction : public FEMaterial
{
public:
  	double	m_ascl;		//!< activation scale factor - lc=".."
	double	m_smax;		//!< max activation
//! force-length curve	
	double	m_ax;
	double	m_ay;
	double	m_bx;
	double	m_by;
	double	m_cx;
	double	m_cy;
	double	m_dx;
	double	m_dy;
	double	m_ex;
	double	m_ey;	
public:
	FENewActiveFiberContraction(FEModel* pfem);
	
	//! initialization
	void Init();

	//! calculate the fiber stress
	double FiberStress(double lamd);

	//! active contraction stiffness contribution
	double FiberStiffness(double lamd);
		
	//! Get the activation
	double GetActivation() { return m_ascl; }	

};

//-----------------------------------------------------------------------------
//! Base class for fiber materials.
class FENewFiberMaterial : public FEMaterial
{
public:
  	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers
	double	m_lam1;		//!< fiber stretch for straightened fibers
public:
	//! Constructor
	FENewFiberMaterial(FEModel* pfem);

	//! Initialization
	void Init();	
		
	//! Calculate the fiber stress
	mat3ds Stress(FEMaterialPoint& mp);

	//! Calculate the fiber tangent
	tens4ds Tangent(FEMaterialPoint& mp);
	
	//! Calculate the fiber strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp);

	void Serialize(DumpFile& ar);
	
	//! Set the active contraction property
	void SetActiveContraction(FENewActiveFiberContraction* pma) { m_pafc = pma; }

	//! get the active contraction property
	FEMaterial* GetActiveContraction() { return m_pafc; }
	
	//! get activation
	double GetActivation() { return (m_pafc ? m_pafc->GetActivation() : 0.0); }
	
	//--- time varying elastance active contraction data ---
	FENewActiveFiberContraction*	m_pafc; // pointer to contraction
	
};