
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "dma.h"
#include "SteffenInterp.h"

/************************************************************************************
 *****	This code was written by Francesca Valsecchi on 07/08/15 to modify MESA WD models.
 *****	The output code is in a format suitable for CAFein.
 *****
 *****	The units used are G = M = R = Rgas = 1
 *****
 *****
 *****	In this new version of the code, I use finite differences to compute derivatives,
 *****	instead of Steffen interpolation
 *****
 *****	http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf
 *****
 *****
 *****
 *****	To set:
 *****	fileNameIn = "log_0.2Msun_cool.data";
 *****	fileNameIn_2 = "logCool.data_0.2Msun.SAIO.data";
 *****	fileNameOut = "MESAToRiccari_log_0.2Msun_cool.data";
 *****	radiusPrecisionInterior = 5.0e-5,
 *****	radiusPrecisionSurface = 1.0e-6,
 *****	MESA_4798 = 1; cause now the header file of profile.data has also the version of MESA, which I can't get rid of.
 *****
 *****	NOTE:
 *****	- I compared the sound speed with the value calculated manually as:
 *****	cSound = sqrt(Gamma1[i]*P[i]/rho[i]) with Gamma calculated by MESA. They agree (1e-8 to 0)
 *************************************************************************************/
using namespace std;

/*************************************************************************************
 *****
 ***** Reading the input model and storing it in a table
 *****
 *************************************************************************************/

int readInputFile(const char *s, vector<vector<double> > &table_func)
{
    ifstream file;
    file.open(s);
    
    int num_lines = 0;
    double data = 0.0;
    while(file)
    {
        string eachLine;
        getline(file, eachLine);
        stringstream elementsOfLine(eachLine);
        eachLine.clear();
        vector<double> eachRow;
        while(elementsOfLine)
        {
            elementsOfLine >> data;
            eachRow.push_back(data);
        }
        if(eachRow.size() > 2) {num_lines++;}
        table_func.push_back(eachRow);
        eachRow.clear();
        elementsOfLine.clear();
    }
    
    
    file.close();
    return num_lines;
}

/*************************************************************************************
 *****
 ***** Manipulation
 *****
 *************************************************************************************/
int
main (void)
{
    cout << setiosflags(ios_base::scientific) << setprecision(16);
    cerr << setiosflags(ios_base::scientific) << setprecision(16);
    
    const int MESA_4798 = 0;
    
    /* the data below for the sun have been taken from http://www.pas.rochester.edu/~emamajek/sun.txt */
    const double PI = 2.0*acos(0.0),
    Rgas = 8.3144621e7,//erg/(K*mol) = erg/(K)
    Gconst_cgs =  6.67259e-8,
    Lsun_cgs = 3.827e33,
    Rsun_cgs  = 6.95660e10,
    Msun_cgs = 1.98855e+33,
    radiusPrecisionInterior = 5.0e-5,
    radiusPrecisionSurface = 5.0e-5;
    
    double Pc = 0.0, Tc = 0.0, rhoc = 0.0, Teff = 0.0, Ltot_cgs = 0.0, Mtot_cgs = 0.0, Rtot_cgs = 0.0,
    rStep = 0.0, dbl_j = 0.0, dummy = 0.0, v_t = 0.0, kappa_ad = 0.0, kappa_s = 0.0,
    eps_ad = 0.0, eps_s = 0.0, dlnLr_dlnr = 0.0,radiusPrecision = 0.0,
    startSurfaceRadiativeRegion = 0.0, endSurfaceRadiativeRegion = 0.0,
    s_units = 0.0, g_units = 0.0, cm_units = 0.0, erg_units = 0.0, K_units = 0.0, U_surf = 0.0;
    
    double *m, *cSound,*Gamma1, *entropy, *r, *rho, *logrho, *T, *logT, *P, *logP, *Lr, *del, *delad, *Cp, *kappa, *g, *N2, *Vg, *U, *V,
    *c1, *c2, *Astar, *L2, *kappa_T, *kappa_rho, *eps_T, *eps_rho, *eps_nucl, *c3, *c4, *scale_heigh, *chiRho, *chiT,
    **fit_m_Coeffs, **fit_rho_Coeffs, **fit_P_Coeffs, **fit_T_Coeffs, **fit_Vg_Coeffs,
    **fit_Astar_Coeffs, **fit_U_Coeffs, **fit_V_Coeffs, **fit_c1_Coeffs, **fit_N2_Coeffs, **fit_L2_Coeffs,
    **fit_g_Coeffs, **fit_Lr_Coeffs, **fit_delad_Coeffs, **fit_del_Coeffs, **fit_Cp_Coeffs, **fit_kappa_Coeffs,
    **fit_kappaRho_Coeffs, **fit_kappaT_Coeffs, **fit_Gamma1_Coeffs, **fit_entropy_Coeffs, **fit_epsT_Coeffs, **fit_epsRho_Coeffs,
    **fit_epsNucl_Coeffs, **fit_c2_Coeffs, **fit_c3_Coeffs , **fit_c4_Coeffs, **fit_scaleHeigh_Coeffs, **fit_chiRho_Coeffs, **fit_chiT_Coeffs,
    *mFuncs, *rhoFuncs, *PFuncs, *TFuncs, *VgFuncs, *AstarFuncs, *UFuncs, *VFuncs, *c1Funcs, *N2Funcs,
    *L2Funcs, *gFuncs, *LrFuncs, *deladFuncs, *delFuncs, *CpFuncs, *kappaFuncs, *kappaRhoFuncs, *kappaTFuncs, *chiRhoFuncs, *chiTFuncs,
    *Gamma1Funcs, *entropyFuncs, *epsTFuncs, *epsRhoFuncs, *epsNuclFuncs, *c2Funcs, *c3Funcs, *c4Funcs, *scale_heighFuncs;
    
    typedef vector<double> Row;
    vector<Row> table, table_2;
    int rowsInput = 0, rowsInput_2 = 0, i = 0, maxMesh = 0, minMesh = 0, j = 0;
    const char *fileNameIn, *fileNameIn_2, *fileNameOut;
    
    /*****************************************************
     *****************************************************
     ** I/O
     *****************************************************
     ****************************************************/
    ifstream input, input_2;
    ofstream output;
    
    fileNameIn = "/Users/francesca/Northwestern/Research/mesaRuns_firstVersionDownloaded/coolMyModel_WD/LOGS/0p23Msun_CAFein/log1673.data";
    fileNameIn_2 = "/Users/francesca/Northwestern/Research/mesaRuns_firstVersionDownloaded/coolMyModel_WD/LOGS/0p23Msun_CAFein/log1673.data.SAIO.txt";
    fileNameOut = "output/WD/0p23Msun/log1673/MESAtoCAFein_log1673_0p23Msun_Z0p02_WD_deltaX5eM5in_07072015.dat";
    
    input.open(fileNameIn);
    input_2.open(fileNameIn_2);
    output.open(fileNameOut);
    output << setiosflags(ios_base::scientific) << setprecision(16);
    /*****************************************************
     *****************************************************
     ** Exctracting first line input model
     **      (global properties)
     *****************************************************
     *****************************************************/
    /*in log#.data the first 3 lines are labels and global properties*/
    string firstLine, labels;
    getline(input, firstLine);
    getline(input, firstLine);
    getline(input, firstLine);
    
    /*store the third line that contains the global properties*/
    istringstream elem_firstLine(firstLine);
    firstLine.clear();
    
    /*loop through the third line and skip the first 6 elements that I don't care about*/
    for (i=1; i<=6+MESA_4798; i++)
    {
        getline(elem_firstLine, labels, ' ');
        elem_firstLine >> dummy;
    }
    
    /*7th, 8th and 9th are useful quantities*/
    getline(elem_firstLine, labels, ' ');
    elem_firstLine >> Teff;
    getline(elem_firstLine, labels, ' ');
    elem_firstLine >> Ltot_cgs;
    getline(elem_firstLine, labels, ' ');
    elem_firstLine >> Rtot_cgs;
    
    /*skip other 8 elements*/
    for (i=1; i<=8; i++)
    {
        getline(elem_firstLine, labels, ' ');
        elem_firstLine >> dummy;
    }
    
    /*store the mass*/
    getline(elem_firstLine, labels, ' ');
    
    elem_firstLine >> Mtot_cgs;
    
//    elem_firstLine.clear();
    
    /*****************************************************
     *****************************************************
     ** Exctracting properties from profile
     ** and storing them in table
     *****************************************************
     *****************************************************/
    
    rowsInput = readInputFile(fileNameIn, table);
    rowsInput = rowsInput - 3; //to cope with the labels!
    cout << "Number rows in input file = " << rowsInput << endl;
    input.close();
    /*****************************************************
     *****************************************************
     ** Exctracting properties from SAIO
     ** and storing them in table_2
     *****************************************************
     *****************************************************/
    rowsInput_2 = readInputFile(fileNameIn_2, table_2);
    cout << "Number rows in SAIO file = " << rowsInput_2 << endl;
    input_2.close();
    
    /*****************************************************
     *****************************************************
     ** allocating vectors
     *****************************************************
     *****************************************************/
    m		= dfun(rowsInput);
    cSound = dfun(rowsInput);
    Gamma1  = dfun(rowsInput);
    entropy = dfun(rowsInput);
    r		= dfun(rowsInput);
    rho		= dfun(rowsInput);
    T		= dfun(rowsInput);
    P		= dfun(rowsInput);
    logrho	= dfun(rowsInput);
    logT	= dfun(rowsInput);
    logP	= dfun(rowsInput);
    Lr		= dfun(rowsInput);
    del		= dfun(rowsInput);
    delad	= dfun(rowsInput);
    Cp		= dfun(rowsInput);
    kappa	= dfun(rowsInput);
    g	    = dfun(rowsInput);
    N2		= dfun(rowsInput);
    L2		= dfun(rowsInput);
    Vg		= dfun(rowsInput);
    U		= dfun(rowsInput);
    V		= dfun(rowsInput);
    c1		= dfun(rowsInput);
    Astar	= dfun(rowsInput);
    kappa_T	= dfun(rowsInput);
    kappa_rho = dfun(rowsInput);
    eps_T	  = dfun(rowsInput);
    eps_rho   = dfun(rowsInput);
    eps_nucl  = dfun(rowsInput);
    c2        = dfun(rowsInput);
    c3        = dfun(rowsInput);
    c4        = dfun(rowsInput);
    scale_heigh   = dfun(rowsInput);
    chiRho	= dfun(rowsInput);
    chiT	= dfun(rowsInput);
    
    
    mFuncs	 = dfun(rowsInput);
    rhoFuncs = dfun(rowsInput);
    PFuncs	 = dfun(rowsInput);
    TFuncs	 = dfun(rowsInput);
    VgFuncs	 = dfun(rowsInput);
    AstarFuncs	= dfun(rowsInput);
    UFuncs	 = dfun(rowsInput);
    VFuncs	 = dfun(rowsInput);
    c1Funcs	 = dfun(rowsInput);
    N2Funcs	 = dfun(rowsInput);
    L2Funcs	 = dfun(rowsInput);
    gFuncs	 = dfun(rowsInput);
    LrFuncs	 = dfun(rowsInput);
    deladFuncs = dfun(rowsInput);
    delFuncs = dfun(rowsInput);
    CpFuncs  = dfun(rowsInput);
    kappaFuncs  = dfun(rowsInput);
    kappaTFuncs    = dfun(rowsInput);
    kappaRhoFuncs  = dfun(rowsInput);
    Gamma1Funcs    = dfun(rowsInput);
    entropyFuncs   = dfun(rowsInput);
    epsTFuncs      = dfun(rowsInput);
    epsRhoFuncs    = dfun(rowsInput);
    epsNuclFuncs   = dfun(rowsInput);
    c2Funcs        = dfun(rowsInput);
    c3Funcs        = dfun(rowsInput);
    c4Funcs        = dfun(rowsInput);
    scale_heighFuncs= dfun(rowsInput);
    chiRhoFuncs  = dfun(rowsInput);
    chiTFuncs    = dfun(rowsInput);
    
    fit_m_Coeffs	 = dfunc(4, rowsInput);
    fit_rho_Coeffs	 = dfunc(4, rowsInput);
    fit_P_Coeffs	 = dfunc(4, rowsInput);
    fit_T_Coeffs	 = dfunc(4, rowsInput);
    fit_Vg_Coeffs	 = dfunc(4, rowsInput);
    fit_Astar_Coeffs = dfunc(4, rowsInput);
    fit_U_Coeffs	 = dfunc(4, rowsInput);
    fit_V_Coeffs	 = dfunc(4, rowsInput);
    fit_c1_Coeffs	 = dfunc(4, rowsInput);
    fit_N2_Coeffs	 = dfunc(4, rowsInput);
    fit_L2_Coeffs	 = dfunc(4, rowsInput);
    fit_g_Coeffs	 = dfunc(4, rowsInput);
    fit_Lr_Coeffs	 = dfunc(4, rowsInput);
    fit_delad_Coeffs = dfunc(4, rowsInput);
    fit_del_Coeffs   = dfunc(4, rowsInput);
    fit_Cp_Coeffs    = dfunc(4, rowsInput);
    fit_kappa_Coeffs     = dfunc(4, rowsInput);
    fit_kappaT_Coeffs    = dfunc(4, rowsInput);
    fit_kappaRho_Coeffs  = dfunc(4, rowsInput);
    fit_Gamma1_Coeffs    = dfunc(4, rowsInput);
    fit_entropy_Coeffs   = dfunc(4, rowsInput);
    fit_epsT_Coeffs      = dfunc(4, rowsInput);
    fit_epsRho_Coeffs    = dfunc(4, rowsInput);
    fit_epsNucl_Coeffs   = dfunc(4, rowsInput);
    fit_c2_Coeffs        = dfunc(4, rowsInput);
    fit_c3_Coeffs        = dfunc(4, rowsInput);
    fit_c4_Coeffs        = dfunc(4, rowsInput);
    fit_scaleHeigh_Coeffs= dfunc(4, rowsInput);
    fit_chiRho_Coeffs     = dfunc(4, rowsInput);
    fit_chiT_Coeffs       = dfunc(4, rowsInput);
    
    /*****************************************************
     *****************************************************
     *****	Reading all the input quantities in cgs units
     *****
     *****	Note, for the SAIO output:
     *****
     *****	[0] i, [1] r, [2] logtau, [3] logT, [4] logrho,
     *****	[5] logP, [6] cv, [7] chiRho, [8] chiT,
     *****	[9] gradT, [10] grada, [11] kap, [12] dlnkap_dlnd,
     *****	[13] dlnkap_dlnT, [14] Lrad, [15] X, [16] Y, [17] L,
     *****	[18] m, [19] eps_nuc, [20] deps_dlnd, [21] deps_dlnT,
     *****	[22] logW
     *****
     *****************************************************
     *****************************************************/
    for(i = 1; i<=rowsInput_2; i++)
    {
        kappa_rho[i-1] = table_2[i][12];//no dimensions
        kappa_T[i-1]   = table_2[i][13];//no dimensions
        
        chiRho[i-1]    = table_2[i][7];//no dimensions
        chiT[i-1]      = table_2[i][8];//no dimensions
        eps_T[i-1]     = table_2[i][21];//no dimensions
        eps_rho[i-1]   = table_2[i][20];//no dimensions
        eps_nucl[i-1]  = table_2[i][19];//erg g-1s-1
    }
    table_2.clear();
    
    for (i=rowsInput+5; i>=6; i--)
    {
        m[rowsInput+5 - i]		= table[i][0]* Msun_cgs; //cgs
        r[rowsInput+5 - i]		= table[i][1]* Rsun_cgs; //cgs
        logrho[rowsInput+5 - i]	= table[i][2];			 //g/cm3
        T[rowsInput+5 - i]		= table[i][3];			 //K
        P[rowsInput+5 - i]		= table[i][4];			 //erg/cm3 = g/(cm s2)
        Lr[rowsInput+5 - i]		= table[i][5]* Lsun_cgs; //cgs
        cSound[rowsInput+5 - i] = table[i][6];			 //cm/s
        del[rowsInput+5 - i]	= table[i][7];			 //no dimensions
        delad[rowsInput+5 - i]	= table[i][8];			 //no dimensions
        Gamma1[rowsInput+5 - i] = table[i][9];			 //no dimensions
        entropy[rowsInput+5 - i] = pow(10.0, table[i][10]);//erg/(K*g) entropy per unit mass
        Cp[rowsInput+5 - i]		= table[i][11];				//Cp in erg/(gK)
        kappa[rowsInput+5 - i]	= table[i][12];				//kappa in cm^2/g
        N2[rowsInput+5 - i]		= table[i][13];				//1/s2
        scale_heigh[rowsInput+5 - i]= table[i][16]* Rsun_cgs;	//in cm
    }
    table.clear();
    
    /*central values from the structure */
    Pc =  P[0];						//cgs
    Tc  = T[0];						//cgs
    rhoc = pow(10.0, logrho[0]);	//cgs
    
    /* Total quantities in cgs units */
    Mtot_cgs = m[rowsInput-1];
    Rtot_cgs = r[rowsInput-1];
    Ltot_cgs = Lr[rowsInput-1];
    
    cout << "**********************************" << endl;
    cout << "******* input model properties" << endl;
    cout << "**********************************\n" << endl;
    cout << "M_cgs ........ = " << Mtot_cgs << endl;
    cout << "M_Msun ........= " << Mtot_cgs/Msun_cgs << endl;
    cout << "R_cgs ........ = " << Rtot_cgs << endl;
    cout << "R_Rsun ....... = " << Rtot_cgs/Rsun_cgs << endl;
    cout << "L_cgs ........ = " << Ltot_cgs << endl;
    cout << "Pc_cgs........ = " << Pc << endl;
    cout << "Tc_cgs........ = " << Tc << endl;
    cout << "rhoc_cgs...... = " << rhoc << endl;
    cout << "log(rhoc_cgs). = " << log10(rhoc) << endl;
    cout << "Teff_csg...... = " << Teff << endl;
    cout << "log(Teff_csg). = " << log10(Teff) << endl;
    cout << "log(L_cgs) ... = " << log10(Ltot_cgs) << endl;
    cout << "\n**********************************" << endl;
    cout << "**********************************" << endl;
    /*****************************************************
     *****************************************************
     ****
     **** Compute the various quantities...
     ****
     *****************************************************
     *****************************************************/
    
    for (i=0; i<rowsInput; i++)
    {
        /**********************************************************************
         At this point all the quantities read are in cgs (or I made them cgs).
         I want to make everything dimensionless.
         Below I define my time, weight, temperature, and length unit. I also
         define the unit of energy for convenience.
         **********************************************************************/
        s_units	  = sqrt(pow(Rtot_cgs, 3.0)/(Gconst_cgs*Mtot_cgs));
        g_units   = Mtot_cgs;
        cm_units  = Rtot_cgs;
        erg_units = Gconst_cgs * pow(Mtot_cgs, 2.0)/Rtot_cgs;
        K_units   = Gconst_cgs * pow(Mtot_cgs, 2.0)/(Rtot_cgs*Rgas);
        /**********************************************************************
         Making everything dimensionless
         **********************************************************************/
        
        /*pressure scale heigh*/
        scale_heigh[i] = scale_heigh[i]/Rtot_cgs;
        
        /* Lamb frequency 1/sec^2*/
        L2[i] = pow(cSound[i]/r[i], 2.0);    //over l(l+1)
        
        /* Lamb frequency dimensionless*/
        L2[i] = L2[i]*pow(s_units,2.0);    //over l(l+1)
        
        /* Radius and mass normalized from 0 to 1*/
        r[i] = r[i]/cm_units;
        m[i] = m[i]/g_units;
        
        /* Gravity */
        g[i] = m[i]/pow(r[i], 2.0);
        
        /* Sound speed */
        cSound[i] = cSound[i]/(cm_units/s_units);
        
        /* Quantity needed by CAFein */
        Vg[i] = r[i]*g[i]/pow(cSound[i], 2.0);
        
        c1[i] = pow(r[i],3.0)/m[i];
        
        /* Specific heat */
        Cp[i] = Cp[i]/(erg_units/(g_units * K_units));
        
        /* Opacity (?)  */
        kappa[i] = kappa[i]/(cm_units*cm_units/g_units);
        
        /* Density */
        logrho[i] = logrho[i] - log10(g_units/pow(cm_units, 3.0));
        rho[i] = pow(10.0, logrho[i]);
        
        /* Temperature */
        T[i] = T[i]/K_units;
        logT[i] = log10(T[i]);
        
        
        /* Pressure */
        P[i] = P[i]/(g_units/(cm_units*pow(s_units, 2.0)));
        logP[i] = log10(P[i]);
        
        
        /* Entropy  */
        entropy[i] = entropy[i]/(erg_units/(g_units * K_units));
        
        
        /* Brunt-Vaisala frequency squared*/
        N2[i] = N2[i]*pow(s_units,2.0);
        
        /* quantity needed by CAFein */
        Astar[i] = r[i]*N2[i]/g[i];
        
        U[i] = 4.0*PI*rho[i]*pow(r[i], 3.0)/m[i];
        
        V[i] = m[i]*rho[i]/(r[i]*P[i]);
        
        
        /*Luminosity*/
        Lr[i] = Lr[i]/(erg_units/s_units);
        
        /* Nuclear energy generation rate */
        eps_nucl[i] = eps_nucl[i]/(erg_units/(g_units*s_units));
        
        
        /* c3 dimensionless */
        c3[i]  = 4.0 * PI * pow(r[i], 3.0)*rho[i]*eps_nucl[i]/Lr[i];
        
        /* c4 dimensionless but calculated with dimensional quantitiess */
        c4[i]  = 4.0 * PI * pow(r[i], 3.0)*rho[i]*T[i]*Cp[i]/Lr[i];
        
    }//for (i=0; i<rowsInput; i++)
    
    U_surf = U[rowsInput-1];
    cout << "surface U  = " << U_surf << endl;
    
    /*finding the size of the radiative region at the very surface (above convective region). This region is present in a hot WD.*/
    for (i=0; i<rowsInput-1; i++)
    {
        if (N2[i] < 0 && N2[i+1] > 0){startSurfaceRadiativeRegion = r[i+1];}
        if (N2[i] > 0 && N2[i+1] < 0){endSurfaceRadiativeRegion = r[i];}
    }
    /***************************************************
     **** Calculate Steffen coefficients...
     **************************************************/
    calculateSteffenInterpolant(rowsInput, r, m, fit_m_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, rho, fit_rho_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, P, fit_P_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, T, fit_T_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, Vg, fit_Vg_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, Astar, fit_Astar_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, U, fit_U_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, V, fit_V_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, c1, fit_c1_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, N2, fit_N2_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, L2, fit_L2_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, g, fit_g_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, Lr, fit_Lr_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, delad, fit_delad_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, del, fit_del_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, Cp, fit_Cp_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, kappa, fit_kappa_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, kappa_rho, fit_kappaRho_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, kappa_T, fit_kappaT_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, Gamma1, fit_Gamma1_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, entropy, fit_entropy_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, eps_T, fit_epsT_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, eps_rho, fit_epsRho_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, eps_nucl, fit_epsNucl_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, c3, fit_c3_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, c4, fit_c4_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, scale_heigh, fit_scaleHeigh_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, chiRho, fit_chiRho_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, chiT, fit_chiT_Coeffs);
    /***************************************************
     **** Compute the missing quantity c2 and calculate its
     **** Steffen interpolants...
     **************************************************/
    
    for (i=0; i<rowsInput; i++)
    {
        kappa_ad = kappa_T[i]*delad[i] + kappa_rho[i]/Gamma1[i];
        
        evaluateSteffenInterpolant(rowsInput,r,fit_delad_Coeffs,r[i], deladFuncs);
        c2[i] = (kappa_ad - 4.0*delad[i])*V[i]*del[i]+delad[i]*((r[i]/delad[i])*deladFuncs[1]+V[i]);
    }
    
    calculateSteffenInterpolant(rowsInput, r, c2, fit_c2_Coeffs);
    
    output << "........M(Msun)..................R(Rsun)..................L(Lsun).................";
    output << "Pc(erg/cm3).................Tc(K).................rhoC(g/cm3).............";
    output << "Teff(K)..................sec_units...............grams_units...............";
    output << "cm_units.................erg_units...............K_units..................unused" << endl;
    
    output << Mtot_cgs/Msun_cgs << "   " << Rtot_cgs/Rsun_cgs << "   " << Ltot_cgs/Lsun_cgs << "   " << Pc << "   " << Tc << "   "
    << rhoc << "   " << Teff << "   "
    << s_units << "   " << g_units << "   " << cm_units << "   " << erg_units << "   " << K_units << "   " << "   0.0000000000000000" << endl;
    
    output << "\n..............m.....................logrho....................logP..........................logT";
    output << ".................r.........................Vg......................Astar.......................";
    output << "U........................c1....................N2....................L2/(l(l+1))....................";
    output << "g.......................Lr.........................V................del_ad(noDim)...............";
    output << "del(noDim)...............Cp...........................v_t.................kappa_s(noDim)..............";
    output << "eps_ad(noDim).........eps_s(noDim)...................c2........................c3.....................";
    output << "c4.....................dlnLr/dlnr...............P_scale..............StartRadSurfLayer..........";
    output << "EndRadSurfLayer.............Gamma1....................entropy.................kappa...................chiRho...................";
    output << "chiT...................U/Usurf..............ddel_ad_dr...............kappa_ad"<< endl;
    
    /* dimensionless v_t*/
    v_t = Cp[0]*delad[0]*rho[0]*T[0]/P[0];
    
    kappa_ad = kappa_T[0]*delad[0] + kappa_rho[0]/Gamma1[0];
    kappa_s = kappa_T[0] - v_t*kappa_rho[0];
    
    evaluateSteffenInterpolant(rowsInput,r,fit_Lr_Coeffs,r[0], LrFuncs);
    dlnLr_dlnr = (r[0]/Lr[0])*LrFuncs[1];
    
    eps_ad = eps_T[0]*delad[0] + eps_rho[0]/Gamma1[0];
    eps_s = eps_T[0] - v_t*eps_rho[0];
    
    evaluateSteffenInterpolant(rowsInput,r,fit_delad_Coeffs,r[0], deladFuncs);
    
    output << "   " << m[0] << "   " << log10(rho[0]) << "   " << log10(P[0]) << "   " <<  log10(T[0]) << "   " << r[0] << "   "
    << Vg[0] << "   " << Astar[0] <<  "   " << U[0] << "   " << c1[0] << "   " << N2[0] << "   "
    << L2[0] << "   " << g[0] << "   " << Lr[0] << "   " << V[0] << "   " << delad[0]<< "   " << del[0] << "   " << Cp[0] << "   "
    << v_t << "   " << kappa_s <<"   " << eps_ad <<"   " << eps_s <<  "   " << c2[0] <<"   " << c3[0] << "   " << c4[0] << "   "
    << dlnLr_dlnr<< "   " << scale_heigh[0] << "   " << startSurfaceRadiativeRegion  << "   " << endSurfaceRadiativeRegion << "   "
    << Gamma1[0] << "   " << entropy[0] << "   " << kappa[0] << "   " << chiRho[0] << "   " << chiT[0] <<  "   " << U[0]/U_surf <<  "   "
    << deladFuncs[1] <<  "   " << kappa_ad << endl;
    
    //Adding data points to see if the fitting with R can help
    const char *out_tempo_address;
    out_tempo_address = "r_delAd_Lr.dat";
    ofstream out_tempo;
    out_tempo.open(out_tempo_address);
    out_tempo << setiosflags(ios_base::scientific) << setprecision(16);
    double constA = 0.0, constB = 0.0;
    
    for (i = 0; i<rowsInput; i++)
    {
        //problematic intervals: 0.74 - 0.85, 0.9996 - 0.99994
        if ((r[i] >= 0.74 && r[i] <= 0.85) || (r[i] >= 0.9996 && r[i] <= 0.99994))
        {
            maxMesh = 3;
            minMesh = 1;
            for (j = minMesh; j <= maxMesh; j++)
            {
                dbl_j = static_cast<double>(j);
                rStep = (dbl_j*r[i] + (maxMesh - dbl_j)*r[i-1])/maxMesh;
                
                evaluateSteffenInterpolant(rowsInput,r,fit_delad_Coeffs,rStep, deladFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_Lr_Coeffs,rStep, LrFuncs);
                
                out_tempo << rStep << "   " << deladFuncs[0] << "   " << deladFuncs[1] << "   " << LrFuncs[0] << "  " << LrFuncs[1] << endl;
            }//for (j = minMesh; j <= maxMesh; j++)
        }//if ((r[i] - r[i-1])>1.0e-4)
        else
        {
            evaluateSteffenInterpolant(rowsInput,r,fit_delad_Coeffs,r[i], deladFuncs);
            evaluateSteffenInterpolant(rowsInput,r,fit_Lr_Coeffs,r[i], LrFuncs);
            
            out_tempo << r[i] << "   " << delad[i] << "   " << deladFuncs[1] << "   " << Lr[i] << "   " << LrFuncs[i] << endl;
        }
    }//for (i = 1; i<rowsInput; i++)
    out_tempo.close();
    
    
    for (i = 1; i<rowsInput; i++)
    {
        /* dimensionless v_t*/
        v_t = Cp[i]*delad[i]*rho[i]*T[i]/P[i];
        
        
        kappa_ad = kappa_T[i]*delad[i] + kappa_rho[i]/Gamma1[i];
        kappa_s = kappa_T[i] - v_t*kappa_rho[i];
        
        evaluateSteffenInterpolant(rowsInput,r,fit_Lr_Coeffs,r[i], LrFuncs);
        dlnLr_dlnr = (r[i]/Lr[i])*LrFuncs[1];
        
        eps_ad = eps_T[i]*delad[i] + eps_rho[i]/Gamma1[i];
        eps_s = eps_T[i] - v_t*eps_rho[i];
        
        if(r[i] < 0.998){radiusPrecision = radiusPrecisionInterior;}
        else {radiusPrecision = radiusPrecisionSurface;}
        
        if ((r[i] - r[i-1])>radiusPrecision)
        {
            maxMesh = static_cast<int>((r[i] - r[i-1])/radiusPrecision);
            minMesh = 1;
            for (j = minMesh; j <= maxMesh; j++)
            {
                dbl_j = static_cast<double>(j);
                
                rStep = (dbl_j*r[i] + (maxMesh - dbl_j)*r[i-1])/maxMesh;
                
                evaluateSteffenInterpolant(rowsInput,r,fit_m_Coeffs,rStep, mFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_rho_Coeffs,rStep, rhoFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_P_Coeffs,rStep, PFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_T_Coeffs,rStep, TFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_Vg_Coeffs,rStep, VgFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_Astar_Coeffs,rStep, AstarFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_U_Coeffs,rStep, UFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_V_Coeffs,rStep, VFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_c1_Coeffs,rStep, c1Funcs);
                evaluateSteffenInterpolant(rowsInput,r,fit_N2_Coeffs,rStep, N2Funcs);
                evaluateSteffenInterpolant(rowsInput,r,fit_L2_Coeffs,rStep, L2Funcs);
                evaluateSteffenInterpolant(rowsInput,r,fit_g_Coeffs,rStep, gFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_Lr_Coeffs,rStep, LrFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_delad_Coeffs,rStep, deladFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_del_Coeffs,rStep, delFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_Cp_Coeffs,rStep, CpFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_kappa_Coeffs,rStep, kappaFuncs);
                
                evaluateSteffenInterpolant(rowsInput,r,fit_kappaRho_Coeffs,rStep, kappaRhoFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_kappaT_Coeffs,rStep, kappaTFuncs);
                
                evaluateSteffenInterpolant(rowsInput,r,fit_Gamma1_Coeffs,rStep, Gamma1Funcs);
                evaluateSteffenInterpolant(rowsInput,r,fit_entropy_Coeffs,rStep, entropyFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_epsT_Coeffs,rStep, epsTFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_epsRho_Coeffs,rStep, epsRhoFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_epsNucl_Coeffs,rStep, epsNuclFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_c2_Coeffs,rStep, c2Funcs);
                evaluateSteffenInterpolant(rowsInput,r,fit_c3_Coeffs,rStep, c3Funcs);
                evaluateSteffenInterpolant(rowsInput,r,fit_c4_Coeffs,rStep, c4Funcs);
                evaluateSteffenInterpolant(rowsInput,r,fit_scaleHeigh_Coeffs,rStep, scale_heighFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_chiRho_Coeffs,rStep, chiRhoFuncs);
                evaluateSteffenInterpolant(rowsInput,r,fit_chiT_Coeffs,rStep, chiTFuncs);
                
                v_t = CpFuncs[0]*deladFuncs[0]*rhoFuncs[0]*TFuncs[0]/PFuncs[0];
                
                kappa_ad = kappaTFuncs[0]*deladFuncs[0] + kappaRhoFuncs[0]/Gamma1Funcs[0];
                kappa_s = kappaTFuncs[0] - v_t*kappaRhoFuncs[0];
                
                dlnLr_dlnr = (rStep/LrFuncs[0])*LrFuncs[1];
                
                eps_ad = epsTFuncs[0]*deladFuncs[0] + epsRhoFuncs[0]/Gamma1Funcs[0];
                eps_s = epsTFuncs[0] - v_t*epsRhoFuncs[0];
                
                evaluateSteffenInterpolant(rowsInput,r,fit_delad_Coeffs,rStep, deladFuncs);
                
                output << "   " << mFuncs[0] << "   " << log10(rhoFuncs[0]) << "   "
                << log10(PFuncs[0]) << "   " <<  log10(TFuncs[0]) << "   " << rStep << "   "
                << VgFuncs[0] << "   " << AstarFuncs[0] <<  "   " << UFuncs[0] << "   "
                << c1Funcs[0] << "   " << N2Funcs[0] << "   " << L2Funcs[0] << "   "
                << gFuncs[0] << "   " << LrFuncs[0] << "   " << VFuncs[0] << "   "
                << deladFuncs[0] << "   " << delFuncs[0] << "   " << CpFuncs[0] <<  "   " << v_t <<  "   "
                << kappa_s<<  "   " << eps_ad <<"   " << eps_s <<  "   " << c2Funcs[0] << "   " << c3Funcs[0] << "   "
                << c4Funcs[0] << "   " << dlnLr_dlnr << "   " << scale_heighFuncs[0]  << "   " << startSurfaceRadiativeRegion  << "   "
                << endSurfaceRadiativeRegion << "   " << Gamma1Funcs[0] << "   " << entropyFuncs[0] << "   " << kappaFuncs[0] << "   "
                << chiRhoFuncs[0] << "   " << chiTFuncs[0] <<  "   " << UFuncs[0]/U_surf <<  "   " << deladFuncs[1] <<  "   " << kappa_ad << endl;
                
            }//for (j = minMesh; j <= maxMesh; j++)
        }//if ((r[i] - r[i-1])>1.0e-4)
        else
        {
            
            evaluateSteffenInterpolant(rowsInput,r,fit_delad_Coeffs,r[i], deladFuncs);
            
            output << "   " << m[i] << "   " << logrho[i] << "   " << logP[i] << "   " <<  logT[i] << "   " << r[i] << "   "
            << Vg[i] << "   " << Astar[i] <<  "   " << U[i] << "   " << c1[i] << "   " << N2[i] << "   "
            << L2[i] << "   " << g[i] << "   " << Lr[i] << "   " << V[i] << "   " << delad[i] << "   " << del[i]<< "   "
            << Cp[i] << "   " << v_t<< "   " << kappa_s<<"   " << eps_ad <<"   " << eps_s <<  "   " << c2[i] << "   "
            << c3[i] << "   " << c4[i]  <<"   " << dlnLr_dlnr << "   " << scale_heigh[i]  << "   " << startSurfaceRadiativeRegion  << "   "
            << endSurfaceRadiativeRegion << "   " << Gamma1[i] << "   " << entropy[i] << "   " << kappa[i] << "   " << chiRho[i] << "   "
            << chiT[i] << "   " << U[i]/U_surf <<  "   " << deladFuncs[1] <<  "   " << kappa_ad << endl;
        }
    }//for (i = 1; i<rowsInput; i++)
    output.close();
    
    
    free(m);
    free(cSound);
    free(Gamma1);
    
    free(r);
    free(rho);
    free(T);
    free(P);
    free(logrho);
    free(logT);
    free(logP);
    free(Lr);
    free(del);
    free(delad);
    free(Cp);
    free(kappa);
    free(g);
    free(N2);
    free(Vg);
    free(U);
    free(V);
    free(c1);
    free(Astar);
    free(L2);
    free(kappa_T);
    free(kappa_rho);
    free(eps_T);
    free(eps_rho);
    free(eps_nucl);
    free(c2);
    free(c3);
    free(c4);
    free(scale_heigh);
    free(chiRho);
    free(chiT);
    
    free(mFuncs);
    free(rhoFuncs);
    free(PFuncs);
    free(TFuncs);
    free(VgFuncs);
    free(AstarFuncs);
    free(UFuncs);
    free(VFuncs);
    free(c1Funcs);
    free(N2Funcs);
    free(L2Funcs);
    free(gFuncs);
    free(LrFuncs);
    free(deladFuncs);
    free(delFuncs);
    free(CpFuncs);
    free(kappaFuncs);
    free(kappaTFuncs);
    free(kappaRhoFuncs);
    free(Gamma1Funcs);
    free(epsTFuncs);
    free(epsRhoFuncs);
    free(epsNuclFuncs);
    free(c2Funcs);
    free(c3Funcs);
    free(c4Funcs);
    free(scale_heighFuncs);
    free(chiRhoFuncs);
    free(chiTFuncs);
    
    for (i=0; i<4; i++)
    {
        free(fit_m_Coeffs[i]);
        free(fit_rho_Coeffs[i]);
        free(fit_P_Coeffs[i]);
        free(fit_T_Coeffs[i]);
        free(fit_Vg_Coeffs[i]);
        free(fit_Astar_Coeffs[i]);
        free(fit_U_Coeffs[i]);
        free(fit_V_Coeffs[i]);
        free(fit_c1_Coeffs[i]);
        free(fit_N2_Coeffs[i]);
        free(fit_L2_Coeffs[i]);
        free(fit_g_Coeffs[i]);
        free(fit_Lr_Coeffs[i]);
        free(fit_delad_Coeffs[i]);
        free(fit_del_Coeffs[i]);
        free(fit_Cp_Coeffs[i]);
        free(fit_kappa_Coeffs[i]);
        free(fit_kappaT_Coeffs[i]);
        free(fit_kappaRho_Coeffs[i]);
        free(fit_Gamma1_Coeffs[i]);
        free(fit_epsT_Coeffs[i]);
        free(fit_epsRho_Coeffs[i]);
        free(fit_epsNucl_Coeffs[i]);
        free(fit_c2_Coeffs[i]);
        free(fit_c3_Coeffs[i]);
        free(fit_c4_Coeffs[i]);
        free(fit_scaleHeigh_Coeffs[i]);
        free(fit_chiRho_Coeffs[i]);
        free(fit_chiT_Coeffs[i]);
    }
    free(fit_scaleHeigh_Coeffs);	
    free(fit_m_Coeffs);
    free(fit_rho_Coeffs);
    free(fit_P_Coeffs);
    free(fit_T_Coeffs);
    free(fit_Vg_Coeffs);
    free(fit_Astar_Coeffs);
    free(fit_U_Coeffs);
    free(fit_V_Coeffs);
    free(fit_c1_Coeffs);
    free(fit_N2_Coeffs);
    free(fit_L2_Coeffs);
    free(fit_g_Coeffs);
    free(fit_Lr_Coeffs);
    free(fit_delad_Coeffs);
    free(fit_del_Coeffs);
    free(fit_Cp_Coeffs);
    free(fit_kappa_Coeffs);
    free(fit_kappaT_Coeffs);
    free(fit_kappaRho_Coeffs);
    free(fit_Gamma1_Coeffs);
    free(fit_epsT_Coeffs);
    free(fit_epsRho_Coeffs);
    free(fit_epsNucl_Coeffs);
    free(fit_c2_Coeffs);
    free(fit_c3_Coeffs);
    free(fit_c4_Coeffs);
    free(fit_chiRho_Coeffs);
    free(fit_chiT_Coeffs);
    
    
    cout << "Done, cazzo" << endl;
    
    return 0;
}


