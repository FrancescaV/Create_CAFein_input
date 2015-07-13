
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
 *****	In this step I create a file where noisy quantities are printed out, to then be used
 *****  in R for smooth interpolation.
 *****
 *****	The units used are G = M = R = Rgas = 1
 *****
 *****
 *****	To set:
 *****	fileNameIn = "log_0.2Msun_cool.data";
 *****	fileNameIn_2 = "logCool.data_0.2Msun.SAIO.data";
 *****	fileNameOut = "MESAToRiccari_log_0.2Msun_cool.data";
 *****	MESA_4798 = 1; cause now the header file of profile.data has also the version of MESA, which I can't get rid of.
 *****
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
    const double Gconst_cgs =  6.67259e-8,
    Lsun_cgs = 3.827e33,
    Rsun_cgs  = 6.95660e10,
    Msun_cgs = 1.98855e+33;
    
    double Teff = 0.0, Ltot_cgs = 0.0, Mtot_cgs = 0.0, Rtot_cgs = 0.0,
    rStep = 0.0, dbl_j = 0.0, dummy = 0.0,
    s_units = 0.0, g_units = 0.0, cm_units = 0.0, erg_units = 0.0;
    
    double *m, *r, *Lr, *delad, **fit_Lr_Coeffs, **fit_delad_Coeffs, *LrFuncs, *deladFuncs;
    
    typedef vector<double> Row;
    vector<Row> table;
    int rowsInput = 0, i = 0, maxMesh = 0, minMesh = 0, j = 0;
    const char *fileNameIn, *out_tempo_address;   
    /*****************************************************
     *****************************************************
     ** I/O
     *****************************************************
     ****************************************************/
    ifstream input;
    ofstream out_tempo;
    
    fileNameIn = "/Users/francesca/Northwestern/Research/mesaRuns_firstVersionDownloaded/coolMyModel_WD/LOGS/0p23Msun_CAFein/log1673.data";
    out_tempo_address = "r_delAd_Lr.dat";
    
    input.open(fileNameIn);
    out_tempo.open(out_tempo_address);
    out_tempo << setiosflags(ios_base::scientific) << setprecision(16);
    /*****************************************************
     *****************************************************
     ** Exctracting first line input model (global properties)
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
     ** allocating vectors
     *****************************************************
     *****************************************************/
    m		= dfun(rowsInput);
    r		= dfun(rowsInput);
    Lr		= dfun(rowsInput);
    delad	= dfun(rowsInput);    
    LrFuncs	   = dfun(rowsInput);
    deladFuncs = dfun(rowsInput);    
    fit_Lr_Coeffs	 = dfunc(4, rowsInput);
    fit_delad_Coeffs = dfunc(4, rowsInput);
    /*****************************************************
     *****************************************************
     *****	Reading input quantities in cgs units
     *****************************************************
     *****************************************************/
    for (i=rowsInput+5; i>=6; i--)
    {
        m[rowsInput+5 - i]		= table[i][0]* Msun_cgs; //cgs
        r[rowsInput+5 - i]		= table[i][1]* Rsun_cgs; //cgs
        Lr[rowsInput+5 - i]		= table[i][5]* Lsun_cgs; //cgs
        delad[rowsInput+5 - i]	= table[i][8];			 //no dimensions
    }
    table.clear();
    
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
    cout << "Teff_csg...... = " << Teff << endl;
    cout << "log(Teff_csg). = " << log10(Teff) << endl;
    cout << "log(L_cgs) ... = " << log10(Ltot_cgs) << endl;
    cout << "\n**********************************" << endl;
    cout << "**********************************" << endl;
    /*****************************************************
     *****************************************************
     **** Compute quantities...
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
        /**********************************************************************
         Making everything dimensionless
         **********************************************************************/
        /* Radius and mass normalized from 0 to 1*/
        r[i] = r[i]/cm_units;
        m[i] = m[i]/g_units;

        /*Luminosity*/
        Lr[i] = Lr[i]/(erg_units/s_units);
    }//for (i=0; i<rowsInput; i++)
    /***************************************************
     **** Calculate Steffen coefficients...
     **************************************************/
    calculateSteffenInterpolant(rowsInput, r, Lr, fit_Lr_Coeffs);
    calculateSteffenInterpolant(rowsInput, r, delad, fit_delad_Coeffs);
    /***************************************************
     **** Looking carefully at stellar model, del_ad 
     **** behaves funny in some points. There, I add more
     **** data points to then fit them with R.
     **** 
     **** In the new output file, I print del_ad and its
     **** derivative, as computed with Steffen.
     **** I'll smooth the fit with R and compare the two     
    **************************************************/
	out_tempo << "              r_over_R			         del_Ad	           	 ";
	out_tempo << "d_del_Ad_dr			           Lr			       d_Lr_dr"<< endl;

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
    
    free(m);
    free(r);
    free(Lr);
    free(delad);
    free(LrFuncs);
    free(deladFuncs);
    
    for (i=0; i<4; i++)
    {
        free(fit_Lr_Coeffs[i]);
        free(fit_delad_Coeffs[i]);
    }
    free(fit_Lr_Coeffs);
    free(fit_delad_Coeffs);
    
    
    cout << "Done, cazzo" << endl;
    
    return 0;
}


