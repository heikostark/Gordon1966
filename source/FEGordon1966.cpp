#include "stdafx.h"
#include "FEGordon1966.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEGordon1966, FEUncoupledMaterial)
	ADD_PARAMETER(m_c1, FE_RANGE_GREATER(0.0), "c1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c2, FE_RANGE_GREATER(0.0), "c2")->setUnits(UNIT_PRESSURE);
	
	ADD_PARAMETER(m_fib.m_c3, FE_RANGE_GREATER(0.0), "c3")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_fib.m_c4, FE_RANGE_GREATER(0.0), "c4")->setUnits(UNIT_NONE);
	ADD_PARAMETER(m_fib.m_c5, FE_RANGE_GREATER(0.0), "c5")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_fib.m_lam1, FE_RANGE_GREATER(0.0), "lam_max")->setUnits(UNIT_NONE); // Maximal length
	
	ADD_PARAMETER(m_fib.m_pafc->m_ascl , FE_RANGE_GREATER_OR_EQUAL(0.0), "ascl" ); // Activation	
	ADD_PARAMETER(m_fib.m_pafc->m_smax , FE_RANGE_GREATER(0.0), "smax" ); // Maximal Force
	
	ADD_PARAMETER(m_fib.m_pafc->m_ax , FE_RANGE_GREATER_OR_EQUAL(0.0), "ax" ); // Force-Length-Curve
	ADD_PARAMETER(m_fib.m_pafc->m_ay , FE_RANGE_GREATER_OR_EQUAL(0.0), "ay" );
	ADD_PARAMETER(m_fib.m_pafc->m_bx , FE_RANGE_GREATER_OR_EQUAL(0.0), "bx" );
	ADD_PARAMETER(m_fib.m_pafc->m_by , FE_RANGE_GREATER_OR_EQUAL(0.0), "by" );
	ADD_PARAMETER(m_fib.m_pafc->m_cx , FE_RANGE_GREATER_OR_EQUAL(0.0), "cx" );
	ADD_PARAMETER(m_fib.m_pafc->m_cy , FE_RANGE_GREATER_OR_EQUAL(0.0), "cy" );
	ADD_PARAMETER(m_fib.m_pafc->m_dx , FE_RANGE_GREATER_OR_EQUAL(0.0), "dx" );
	ADD_PARAMETER(m_fib.m_pafc->m_dy , FE_RANGE_GREATER_OR_EQUAL(0.0), "dy" );
	ADD_PARAMETER(m_fib.m_pafc->m_ex , FE_RANGE_GREATER_OR_EQUAL(0.0), "ex" );
	ADD_PARAMETER(m_fib.m_pafc->m_ey , FE_RANGE_GREATER_OR_EQUAL(0.0), "ey" );
	
	ADD_PROPERTY(m_fiber, "fiber")->SetDefaultType("vector");

//	ADD_PROPERTY(m_ac, "active_contraction", FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEGordon1966::FEGordon1966(FEModel* pfem) : FEUncoupledMaterial(pfem), m_fib(pfem)  
{   
	printf("Gordon1966"); 
	m_c1 = 0.0;
	m_c2 = 0.0;

	//m_ac = nullptr;
	//m_fib.SetParent(this);
	m_fiber = nullptr;
}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPointData* FEGordon1966::CreateMaterialPointData() 
{
    // create the elastic solid material point
    FEMaterialPointData* ep = new FEElasticMaterialPoint;
    
    // create the material point from the active contraction material
    //if (m_ac) ep->Append(m_ac->CreateMaterialPointData());

	return ep;
}

//-----------------------------------------------------------------------------
// Data initialization
bool FEGordon1966::Init()
{
	if (FEUncoupledMaterial::Init() == false) return false;
	if (m_fib.Init() == false) return false;	// first init for m_fib.m_pafc
		//if (m_E <= 0) throw MaterialError("Invalid value for E");
	//if (!IN_RIGHT_OPEN_RANGE(m_v, -1.0, 0.5)) throw MaterialRangeError("v", -1.0, 0.5, true, false);
	printf("Muscle force l=%f (%E * %E * (%E,%E,%E,%E,%E))\n",m_fib.m_lam1,m_fib.m_pafc->m_ascl,m_fib.m_pafc->m_smax,m_c1,m_c2,m_fib.m_c3,m_fib.m_c4,m_fib.m_c5);

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEGordon1966::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// material axes
	mat3d Q = GetLocalCS(mp);
	// get the fiber vector in local coordinates
	vec3d fiber = m_fiber->unitVector(mp);
	// convert to global coordinates
	vec3d a0 = Q * fiber;	
	
	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1 = m_c1;
	double W2 = m_c2;
	// ------------------------------------------------

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	// calculate stress s = pI + 2/J * dev(T) 
	mat3ds s = T.dev()*(2.0/J);

	// add the fiber stress
	s += m_fib.Stress(mp,a0);

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEGordon1966::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0/J;

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// material axes
	mat3d Q = GetLocalCS(mp);
	// get the fiber vector in local coordinates
	vec3d fiber = m_fiber->unitVector(mp);
	// convert to global coordinates
	vec3d a0 = Q * fiber;	
	
	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1, W2;
	W1 = m_c1;
	W2 = m_c2;
	// ------------------------------------

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds devs = pt.m_s.dev();

	// mean pressure
	//double p = pt.avgp;

	// d2W/dCdC:C
	mat3ds WCCxC = B*(W2*I1) - B2*W2;

	tens4ds cw = (BxB - B4)*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);
	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	return c + m_fib.Tangent(mp,a0);
}
