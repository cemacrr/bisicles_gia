#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

/// Calving model based on viscous stress 
/// see (for 1D) Benn (2007), Earth-Science Reviews vol 82 p 143 doi:10.1016/j.earscirev.2007.02.002
/// Nick et al (2011) J. Glaciol vol 56 p 781 doi:10.3189/002214310794457344
/// author : Andrew Taylor, University of Bristol

#include "CalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "BennCalvingModel.H"
#include "NamespaceHeader.H"
//#define XX 0
//#define YX 1
//#define XY 2
//#define YY 3



void BennCalvingModel::endTimeStepModifyState
(LevelData<FArrayBox>& a_thickness, 
 const AmrIce& a_amrIce,
 int a_level)
{

  //Real time = a_amrIce.time();
  const Real rhoi = 918.0;
  //const Real rhoo = 1028.0;
  const Real grav = 9.81;

  /// Get viscous tensor (not done yet!)....
  /// LevelData<FArrayBox>& viscousTensor = a_amrIce.viscousTensor();
  //bool calvingActive = (time >= m_startTime && time < m_endTime);
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);



  const LevelData<FArrayBox>& vt  = *a_amrIce.viscousTensor(a_level);



  for (DataIterator dit (levelCoords.grids()); dit.ok(); ++dit)
    {
      const FArrayBox& VT = vt[dit];
      FArrayBox& thck = a_thickness[dit];
      const FArrayBox& topg = levelCoords.getTopography()[dit];
      const FArrayBox& usrf = levelCoords.getSurfaceHeight()[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      //const FArrayBox& usrf = std::max (topg+thck, thck*(1.0-(rhoi/rhoo)))[dit];
      
      // pull box out with no ghost data using [dit] and levelCoords
      const Box& b = levelCoords.grids()[dit];
      //const Box& b = thck.box();
  
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  const int& plus_i = iv[0]+1;
	  const int& minus_i = iv[0]-1;
	  const int& plus_j = iv[1]+1;
	  const int& minus_j = iv[1]-1;
	  const IntVect iv_plus_i(plus_i, iv[1]);
	  const IntVect iv_minus_i (minus_i, iv[1]);
	  const IntVect iv_plus_i_plus_j (plus_i,plus_j);
	  const IntVect iv_plus_i_minus_j (plus_i,minus_j);
	  const IntVect iv_minus_i_plus_j (minus_i,plus_j);
	  const IntVect iv_minus_i_minus_j (minus_i,minus_j);
	  const IntVect iv_plus_j (iv[0], plus_j);
	  const IntVect iv_minus_j (iv[0], minus_j);

	  const int& plus_2i = iv[0]+2;
	  const int& minus_2i = iv[0]-2;
	  const int& plus_2j = iv[1]+2;
	  const int& minus_2j = iv[1]-2;
	  const IntVect iv_plus_2i (plus_2i, iv[1]);
	  const IntVect iv_minus_2i (minus_2i, iv[1]);
	  const IntVect iv_plus_2j (iv[0], plus_2j);
	  const IntVect iv_minus_2j(iv[0], minus_2j);
	  const IntVect iv_plus_2i_plus_2j(plus_2i,plus_2j);
	  const IntVect iv_plus_2i_minus_2j (plus_2i,minus_2j);
	  const IntVect iv_minus_2i_plus_2j (minus_2i,plus_2j);
	  const IntVect iv_minus_2i_minus_2j (minus_2i,minus_2j);

	  

	  const Real& vt_xx = VT(iv,0);

	  const Real& vt_xy = VT(iv,2);
	  const Real& vt_yy = VT(iv,3);
	  const Real& vt_yx = VT(iv,1);

	  
	  // convert VT to Pa if thck > 0.0 and calculate the First Principal
	  // Stress
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck (iv) = 0;
	    }
	  else if (thck(iv) > 0.0)
	    {
	      //Real xx_Pa = vt_xx/thck(iv);
	      //Real xy_Pa = vt_xy/thck(iv);
	      //Real yy_Pa = vt_yy/thck(iv);
	      //Real yx_Pa = vt_yx/thck(iv);
	      //pout() << "converted VT comps to Pa" << std::endl;
	      //Real xy_Pa(iv) = VT(iv, XY)/thck[bit](iv);
	      //Real yy_Pa(iv) = VT(iv, YY)/thck(iv);
	      //Real yx_Pa(iv) = VT(iv, YX)/thck(iv);
	      Real FPS = (((vt_xx + vt_yy)/2) + std::sqrt((((vt_xx-vt_yy)*(vt_xx-vt_yy)/2))+(((vt_xy+vt_yx)/2)*(vt_xy+vt_yx)/2)))/thck(iv);
	      // Calculate crevasse penetration depth using Ds = R / grav * rhoi
	      Real Ds = FPS / (grav*rhoi);
	      
	      // Real usrf = std::max (topg(iv)+thck(iv), thck(iv)*1.0-(rhoi/rhoo));

	      bool open_sea_plus_i = mask(iv_plus_i) == OPENSEAMASKVAL;
	      bool open_sea_plus_j = mask(iv_plus_j) == OPENSEAMASKVAL;
	      bool open_sea_minus_i = mask(iv_minus_i) == OPENSEAMASKVAL;
	      bool open_sea_minus_j = mask(iv_minus_j) == OPENSEAMASKVAL;
	      bool open_sea_adj = open_sea_plus_i || open_sea_plus_j || open_sea_minus_i || open_sea_minus_j;

	      bool open_sea_plus_2i = mask(iv_plus_2i) == OPENSEAMASKVAL;
	      bool open_sea_plus_2j = mask(iv_plus_2j) == OPENSEAMASKVAL;
	      bool open_sea_minus_2i = mask(iv_minus_2i) == OPENSEAMASKVAL;
	      bool open_sea_minus_2j = mask(iv_minus_2j) == OPENSEAMASKVAL;
	      bool open_sea_plus_2i_plus_2j = mask(iv_plus_2i_plus_2j) == OPENSEAMASKVAL;
	      bool open_sea_plus_2i_minus_2j = mask(iv_plus_2i_minus_2j) == OPENSEAMASKVAL;
	      bool open_sea_minus_2i_plus_2j = mask(iv_minus_2i_plus_2j) == OPENSEAMASKVAL;
	      bool open_sea_minus_2i_minus_2j = mask(iv_minus_2i_minus_2j) == OPENSEAMASKVAL;

	      // ALL +2 open sea
	      // Isolated island case
	      bool open_sea_all_adj_isol = open_sea_plus_2i && open_sea_plus_2j && open_sea_minus_2i && open_sea_minus_2j;
	      // Attached case
	      bool open_sea_2i_plus_2i_plus_2j_att = open_sea_plus_2i && open_sea_plus_2j && open_sea_minus_2i &&
		open_sea_minus_2j && open_sea_plus_2i_minus_2j && open_sea_minus_2i_plus_2j && 
		open_sea_minus_2i_minus_2j;
	      bool open_sea_2i_plus_2i_minus_2j_att = open_sea_plus_2i && open_sea_plus_2j && open_sea_minus_2i &&
		open_sea_minus_2j && open_sea_plus_2i_plus_2j && open_sea_minus_2i_plus_2j && 
		open_sea_minus_2i_minus_2j;
	      bool open_sea_2i_minus_2i_plus_2j_att = open_sea_plus_2i && open_sea_plus_2j && open_sea_minus_2i &&
		open_sea_minus_2j && open_sea_plus_2i_plus_2j && open_sea_plus_2i_minus_2j  && 
		open_sea_minus_2i_minus_2j;
	      bool open_sea_2i_minus_2i_minus_2j_att = open_sea_plus_2i && open_sea_plus_2j && open_sea_minus_2i &&
		open_sea_minus_2j && open_sea_plus_2i_plus_2j && open_sea_plus_2i_minus_2j && open_sea_minus_2i_plus_2j;
	      // THREE +2 open sea
	      bool open_sea_three_v1_adj = open_sea_plus_2i && open_sea_plus_2j && open_sea_minus_2i;
	      bool open_sea_three_v2_adj = open_sea_plus_2j && open_sea_minus_2i && open_sea_minus_2j;
	      bool open_sea_three_v3_adj = open_sea_plus_2i && open_sea_minus_2i && open_sea_minus_2j;
	      bool open_sea_three_v4_adj = open_sea_plus_2i && open_sea_plus_2j && open_sea_minus_2j;
	      // TWO CORNER +2 open sea
	      bool open_sea_two_v1_adj = open_sea_plus_2i && open_sea_plus_2j;
	      bool open_sea_two_v2_adj = open_sea_plus_2i && open_sea_minus_2j;
	      bool open_sea_two_v3_adj = open_sea_plus_2j && open_sea_minus_2i;
	      bool open_sea_two_v4_adj = open_sea_minus_2j && open_sea_minus_2i;
	      // ONE +2 open sea
	      bool open_sea_one_plus_2i_adj = open_sea_plus_2i;
	      bool open_sea_one_plus_2j_adj = open_sea_plus_2j;
	      bool open_sea_one_minus_2i_adj = open_sea_minus_2i;
	      bool open_sea_one_minus_2j_adj = open_sea_minus_2j;
	      

	      //bool open_sea_plus_2j_adj = open_sea_plus_2j 
	      //bool open_sea_minus_2i_adj = open_sea_minus_2i;
	      //bool open_sea_minus_2j_adj = open_sea_minus_2j;

	      if ((usrf(iv) - Ds < -0.01) && open_sea_adj) 
		{
		  thck(iv) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_all_adj_isol)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		  thck(iv_plus_i_plus_j) = 0;
		  thck(iv_minus_i_plus_j) = 0;
		  thck(iv_plus_i_minus_j) = 0;
		  thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_2i_plus_2i_plus_2j_att)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		  //thck(iv_plus_i_plus_j) = 0;
		  thck(iv_minus_i_plus_j) = 0;
		  thck(iv_plus_i_minus_j) = 0;
		  thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_2i_plus_2i_minus_2j_att)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		  thck(iv_plus_i_plus_j) = 0;
		  thck(iv_minus_i_plus_j) = 0;
		  //thck(iv_plus_i_minus_j) = 0;
		  thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_2i_minus_2i_plus_2j_att)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		  thck(iv_plus_i_plus_j) = 0;
		  //thck(iv_minus_i_plus_j) = 0;
		  thck(iv_plus_i_minus_j) = 0;
		  thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_2i_minus_2i_minus_2j_att)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		  thck(iv_plus_i_plus_j) = 0;
		  thck(iv_minus_i_plus_j) = 0;
		  thck(iv_plus_i_minus_j) = 0;
		  //thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_three_v1_adj)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  //thck(iv_minus_j) = 0;
		  thck(iv_plus_i_plus_j) = 0;
		  thck(iv_minus_i_plus_j) = 0;
		  //thck(iv_plus_i_minus_j) = 0;
		  //thck(iv_minus_i_minus_j) = 0;
		}

	      else if ((usrf(iv) - Ds < -0.01) && open_sea_three_v2_adj)
		{
		  thck(iv) = 0;
		  //thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		  //thck(iv_plus_i_plus_j) = 0;
		  thck(iv_minus_i_plus_j) = 0;
		  //thck(iv_plus_i_minus_j) = 0;
		  thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_three_v3_adj)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  //thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		  //thck(iv_plus_i_plus_j) = 0;
		  //thck(iv_minus_i_plus_j) = 0;
		  thck(iv_plus_i_minus_j) = 0;
		  thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_three_v4_adj)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  //thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		  thck(iv_plus_i_plus_j) = 0;
		  //thck(iv_minus_i_plus_j) = 0;
		  thck(iv_plus_i_minus_j) = 0;
		  //thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_two_v1_adj)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  //thck(iv_minus_i) = 0;
		  //thck(iv_minus_j) = 0;
		  thck(iv_plus_i_plus_j) = 0;
		  //thck(iv_minus_i_plus_j) = 0;
		  //thck(iv_plus_i_minus_j) = 0;
		  //thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_two_v2_adj)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  //thck(iv_plus_j) = 0;
		  //thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		  //thck(iv_plus_i_plus_j) = 0;
		  //thck(iv_minus_i_plus_j) = 0;
		  thck(iv_plus_i_minus_j) = 0;
		  //thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_two_v3_adj)
		{
		  thck(iv) = 0;
		  //thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  //thck(iv_minus_j) = 0;
		  //thck(iv_plus_i_plus_j) = 0;
		  thck(iv_minus_i_plus_j) = 0;
		  //thck(iv_plus_i_minus_j) = 0;
		  //thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_two_v4_adj)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  //thck(iv_minus_i) = 0;
		  //thck(iv_minus_j) = 0;
		  thck(iv_plus_i_plus_j) = 0;
		  thck(iv_minus_i_plus_j) = 0;
		  thck(iv_plus_i_minus_j) = 0;
		  thck(iv_minus_i_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_one_plus_2i_adj)
		{
		  thck(iv) = 0;
		  thck(iv_plus_i) = 0;
		  //thck(iv_plus_j) = 0;
		  //thck(iv_minus_i) = 0;
		  //thck(iv_minus_j) = 0;
		  
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_one_plus_2j_adj)
		{
		  thck(iv) = 0;
		  //thck(iv_plus_i) = 0;
		  thck(iv_plus_j) = 0;
		  //thck(iv_minus_i) = 0;
		  //thck(iv_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_one_minus_2i_adj)
		{
		  thck(iv) = 0;
		  //thck(iv_plus_i) = 0;
		  //thck(iv_plus_j) = 0;
		  thck(iv_minus_i) = 0;
		  //thck(iv_minus_j) = 0;
		}
	      else if ((usrf(iv) - Ds < -0.01) && open_sea_one_minus_2j_adj)
		{
		  thck(iv) = 0;
		  //thck(iv_plus_i) = 0;
		  //thck(iv_plus_j) = 0;
		  //thck(iv_minus_i) = 0;
		  thck(iv_minus_j) = 0;
		} 
	    }
	}
    }
}

#include "NamespaceFooter.H"	      
	      
	  

	  
