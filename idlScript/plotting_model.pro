PRO plotting_model
;poster
;  !P.THICK = 6.8
;  !P.CHARTHICK = 6.8
;  !X.THICK = 6.8
;  !Y.THICK = 6.8
;paper  
 !P.THICK = 3.0
  !P.CHARTHICK = 4.0
  !X.THICK = 3.0
  !Y.THICK = 3.0
; Define default character size
!P.CHARSIZE = 1.0
; Arrange the P.MULTI system varable so that we get four plots per page
; !P.MULTI=[4,2,2]
!X.MARGIN=[5,1]
!Y.MARGIN=[2,2]
;Plot positions
plot1 = [0.11,0.12,0.975,0.96]

MsunToG = 1.98892e33 
RsunToCm = 6.96e10
pi = 3.14159265
thisLetter = "156B
ldeg = 2.0
Compile_Opt DEFINT32


READCOL, "/Users/francesca/Northwestern/Research/create_models_for_tides_and_oscillations/createWDmodelsForOscillationCode_FromMESA_fixedUnitsIn_04022013/output/WD/0p23Msun/log1673/MESAtoCAFein_log1673_0p23Msun_Z0p02_WD_deltaX5eM5in_09222013_toPlot.dat",$
m, logrho, logP, logT, csi, Vg, Astar, U, c1, BVfreq,$ 
lamb_freq, g, Lr, V, del_ad, del, Cp, v_t, kappa_s, eps_ad,$ 
eps_s, c2, c3, c4, dlnLr_dlnr, P_scale, StartRadSurfLayer, EndRadSurfLayer, Gamma1, entropy, $
kappa, chiRho, chiT, UtoUsurf, $
FORMAT = 'D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D'


PS_OPEN, "/Users/francesca/Northwestern/Research/create_models_for_tides_and_oscillations/createWDmodelsForOscillationCode_FromMESA_fixedUnitsIn_04022013/output/WD/0p23Msun/log1673/plots/MESAtoCAFein_log1673_0p23Msun_Z0p02_WD_deltaX5eM5in_09222013_Vg_Astar_U_c1_rho_T.eps",$
          RATIO =1.4, /HICOLOUR
;; ;Load color table
red   = [0,255,195,135, 100, 255, 140, 255,   0,   0,   0,   0,   0,    0, 255, 140, 255, 150]
green = [0,255,195,135, 100, 0  ,   0, 150, 255, 130,   0,   0, 255,  140,   0,   0, 255, 150]
blue  = [0,255,195,135, 100, 0  ,   0,   0,   0,   0, 255, 130, 255,  140, 255, 140,   0,   0]
TVLCT, [red], [green], [blue]
; 0:black, 1:white, 2:light-grey, 3:grey, 4:dark-grey, 5:red, 6:dark- red, 7:orange, 8:green, 9:dark-green
; 10:blue, 11:dark-blue, 12:cyan, 13:dark-cyan, 14:magenta,  15:purple, 16:yellow, 17:dark-yellow

xminPlotAll = 0.1
ymaxPlotAll = 0.99

deltaXplot = 0.45
deltaYplot = 0.13
xCut = 0.99
xmaxCut = 1.000


xminPlot = xminPlotAll
ymaxPlot = ymaxPlotAll

yminPlot = ymaxPlot-deltaYplot
xmaxPlot = xminPlot+deltaXplot

xmin = min(csi)-0.1
xmax = max(csi)+0.1

xmin = 0
xmax = xCut

ymin = min(Vg)-10
ymax = max(Vg)+10


yminVg = 1e-4
ymaxVg = 1e8

PLOT, csi, Vg, /ylog, yRange = [yminVg,ymaxVg],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = 'V!Dg!N',xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, Vg_noMesh, linestyle=1, color=5

ymin = min(Astar)
ymax = max(Astar)
yminAstar = 1e-6
ymaxAstar = 1e6
print,  "A = ", ymin, ymax
ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, Astar, /ylog, yRange = [yminAstar,ymaxAstar],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = 'A!U*!N',xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, Astar_noMesh, linestyle=1, color=5
yminU = min(U)-10
ymaxU = max(U)+10

;yminU = 1e-11
;ymaxU = 1e3
ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, U, /ylog, yRange = [yminU,ymaxU],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = "U",xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, U_noMesh, linestyle=1, color=5


ymin = min(c1)-10
ymax = max(c1)+10
ymaxc1 = 10
yminc1 = 1e-4
ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, c1, /ylog, yRange = [yminc1,ymaxc1],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1,Ytitle = "c!D1!N",xtickformat = '(A1)', $
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, c1_noMesh, linestyle=1, color=5

yminT = min(logT)-1
ymaxT = max(logT)+1

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, logT, yRange = [yminT,ymaxT],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = 'log(T)',xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, logT_noMesh, linestyle=1, color=5


yminRho = min(logRho)-1
ymaxRho = max(logRho)+1

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, logRho, yRange = [yminRho,ymaxRho],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = 'log(!7q!3)',xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, logRho_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

yminc3 = min(c3)
ymaxc3 = max(c3)

yminc3 = -0.5
ymaxc3 = 0.5
PLOT, csi, c3, /ylog, yRange = [yminc3,ymaxc3],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = 'c3', Xtitle = "r", $
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, logRho_noMesh, linestyle=1, color=5


;&&&&&&&&&&&&&&&&&&&&&&&& surface region
deltaXplot = 0.45
xminPlot = xmaxPlot
ymaxPlot = ymaxPlotAll

yminPlot = ymaxPlot-deltaYplot
xmaxPlot = xminPlot+deltaXplot

xmin = xmax
xmax = xmaxCut

PLOT, csi, Vg, /ylog, yRange = [yminVg,ymaxVg],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1,xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, Vg_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, Astar, /ylog, yRange = [yminAstar,ymaxAstar],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, Astar_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, U, /ylog, yRange = [yminU,ymaxU],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, U_noMesh, linestyle=1, color=5


ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, c1, /ylog, yRange = [yminc1,ymaxc1],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1,ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase, xtickformat = '(A1)'
;OPLOT, csi_noMesh, c1_noMesh, linestyle=1, color=5


ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, logT, yRange = [yminT,ymaxT],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1,xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, logT_noMesh, linestyle=1, color=5


ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, logRho, yRange = [yminRho,ymaxRho],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, ytickformat = '(A1)',xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, logRho_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, c3, /ylog, yRange = [yminc3,ymaxc3],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Xtitle = "r", ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, logRho_noMesh, linestyle=1, color=5

ps_close 



PS_OPEN, "/Users/francesca/Northwestern/Research/create_models_for_tides_and_oscillations/createWDmodelsForOscillationCode_FromMESA_fixedUnitsIn_04022013/output/WD/0p23Msun/log1673/plots/MESAtoCAFein_log1673_0p23Msun_Z0p02_WD_deltaX5eM5in_09222013_nonAdiabatic_zoomSurface.eps",$
          RATIO =1.4, /HICOLOUR
;; ;Load color table
red   = [0,255,195,135, 100, 255, 140, 255,   0,   0,   0,   0,   0,    0, 255, 140, 255, 150]
green = [0,255,195,135, 100, 0  ,   0, 150, 255, 130,   0,   0, 255,  140,   0,   0, 255, 150]
blue  = [0,255,195,135, 100, 0  ,   0,   0,   0,   0, 255, 130, 255,  140, 255, 140,   0,   0]
TVLCT, [red], [green], [blue]
; 0:black, 1:white, 2:light-grey, 3:grey, 4:dark-grey, 5:red, 6:dark- red, 7:orange, 8:green, 9:dark-green
; 10:blue, 11:dark-blue, 12:cyan, 13:dark-cyan, 14:magenta,  15:purple, 16:yellow, 17:dark-yellow


;xminPlotAll = 0.1
;ymaxPlotAll = 0.99
xminPlot = xminPlotAll
ymaxPlot = ymaxPlotAll
deltaXplot = 0.45
deltaYplot = 0.12

yminPlot = ymaxPlot-deltaYplot
xmaxPlot = xminPlot+deltaXplot

xmin = 0
xmax = xCut

ymin_c2 = -5e3
ymax_c2 = 1e4

PLOT, csi, c2, yRange = [ymin_c2,ymax_c2],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = 'c!D2!N',xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, c2_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

yminV = 1e-3
ymaxV = 5e4
PLOT, csi, V, /ylog, yRange = [yminV,ymaxV],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = 'V',xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, V_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

yminDelAd = 0.03
ymaxDelAd = 1  

PLOT, csi, del_ad, /ylog, yRange = [yminDelAd,ymaxDelAd],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = '!7D!3(-), !7D!3!Dad!N(...)',xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
oplot, csi, del, linestyle = 1
;OPLOT, csi_noMesh, del_ad_noMesh, linestyle=2, color=5
;OPLOT, csi_noMesh, del_noMesh, linestyle=2, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

yminVt = 0
ymaxVt = 3
PLOT, csi, v_t, yRange = [yminVt,ymaxVt],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = "v!Dt!N",xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, v_t_noMesh, linestyle=1, color=5


ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

yminKs = min(kappa_s)-5
ymaxKs = max(kappa_s)+5
PLOT, csi, kappa_s, yRange = [yminKs,ymaxKs],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = "k!D!Ns",xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, kappa_s_noMesh, linestyle=1, color=5


ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

ymin_eps_ad = min(eps_ad)
ymax_eps_ad = max(eps_ad)
ymin_eps_ad = -10
ymax_eps_ad = 20
PLOT, csi, eps_ad, yRange = [ymin_eps_ad,ymax_eps_ad],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = "!7e!3!Dad!N",xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, eps_ad_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot


ymin_eps_s = -1e-5
ymax_eps_s = 1e-5
ymin_eps_s = -50
ymax_eps_s = 200

PLOT, csi, eps_s, yRange = [ymin_eps_s,ymax_eps_s],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = "!7e!3!Ds!N",xtickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, eps_s_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

ymaxDLr = 1e-3
yminDlr = -1e-3
ymaxDLr = -3
yminDlr = 5
PLOT, csi, dlnLr_dlnr, yRange = [yminDLr,ymaxDLr],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1,Xtitle = "!4"+String(thisLetter)+"!X", Ytitle = "dln Lr/dln r", $
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, dlnLr_dlnr_noMesh, linestyle=1, color=5

;&&&&&&&&&&&&&&&&&&&&&&&& surface region
deltaXplot = 0.45
xminPlot = xmaxPlot
ymaxPlot = ymaxPlotAll

yminPlot = ymaxPlot-deltaYplot
xmaxPlot = xminPlot+deltaXplot

xmin = xmax
xmax = xmaxCut


PLOT, csi, c2, yRange = [ymin_c2,ymax_c2],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, c2_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, V, /ylog, yRange = [yminV,ymaxV],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, V_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, del_ad, /ylog, yRange = [yminDelAd,ymaxDelAd],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, del_ad_noMesh, linestyle=2, color=5
;OPLOT, csi_noMesh, del_noMesh, linestyle=2, color=5
oplot, csi, del, linestyle = 1

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot


PLOT, csi, v_t, yRange = [yminVt,ymaxVt],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, v_t_noMesh, linestyle=1, color=5


ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, kappa_s, yRange = [yminKs,ymaxKs],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, kappa_s_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, eps_ad, yRange = [ymin_eps_ad,ymax_eps_ad],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, eps_ad_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, eps_s, yRange = [ymin_eps_s,ymax_eps_s],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, xtickformat = '(A1)',ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, eps_s_noMesh, linestyle=1, color=5

ymaxPlot = yminPlot
yminPlot = ymaxPlot-deltaYplot

PLOT, csi, dlnLr_dlnr, yRange = [yminDlr,ymaxDlr],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1,Xtitle = "r",ytickformat = '(A1)',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
;OPLOT, csi_noMesh, dlnLr_dlnr_noMesh, linestyle=1, color=5
ps_close 

PS_OPEN, "/Users/francesca/Northwestern/Research/create_models_for_tides_and_oscillations/createWDmodelsForOscillationCode_FromMESA_fixedUnitsIn_04022013/output/WD/0p23Msun/log1673/plots/MESAtoCAFein_log1673_0p23Msun_Z0p02_WD_deltaX5eM5in_09222013_BV_lamb_zoomSurface.eps",$
          RATIO =1.0, /HICOLOUR
;; ;Load color table
red   = [0,255,195,135, 100, 255, 140, 255,   0,   0,   0,   0,   0,    0, 255, 140, 255, 150]
green = [0,255,195,135, 100, 0  ,   0, 150, 255, 130,   0,   0, 255,  140,   0,   0, 255, 150]
blue  = [0,255,195,135, 100, 0  ,   0,   0,   0,   0, 255, 130, 255,  140, 255, 140,   0,   0]
TVLCT, [red], [green], [blue]
; 0:black, 1:white, 2:light-grey, 3:grey, 4:dark-grey, 5:red, 6:dark- red, 7:orange, 8:green, 9:dark-green
; 10:blue, 11:dark-blue, 12:cyan, 13:dark-cyan, 14:magenta,  15:purple, 16:yellow, 17:dark-yellow

xminPlot = 0.1
xmaxPlot = 0.9
ymaxPlot = 0.95
deltaYplot = 0.4
yminPlot = ymaxPlot-deltaYplot

xmin = -0.1
xmax = 1.1
ymin = min(BVfreq)
ymax = max(BVfreq)
print, "BV extremes?", ymin, ymax
ymin = 1e-4
ymax = 1e3
PLOT, csi, BVfreq,/ylog,yRange = [ymin,ymax],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1,Ytitle = 'N!U2!N and L!D2!N!U2!N/GMR!U-3!N',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
OPLOT, csi,lamb_freq*ldeg*(ldeg+1.0), linestyle = 1;
;OPLOT, csi_noMesh, BVfreq_noMesh, linestyle=1, color=5
;OPLOT, csi_noMesh, lamb_freq_noMesh*ldeg*(ldeg+1.0), linestyle=1, color=5



ymaxPlot = yminPlot-deltaYplot/5
yminPlot = ymaxPlot-deltaYplot

ymin = 1e-3
ymax = 1e1
;xmin = xCut
;xmax = XmaxCut
PLOT, csi, BVfreq,/ylog,yRange = [ymin,ymax],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1,Xtitle = '!4'+String(thisLetter)+'!X/!4'+String(thisLetter)+'!X!D1!N', Ytitle = 'N!U2!N and L!D2!N!U2!N/GMR!U-3!N',$
position = [xminPlot, yminPlot, xmaxPlot, ymaxPlot], /noerase
OPLOT, csi,lamb_freq*ldeg*(ldeg+1.0), linestyle = 1;
;plots, [xmin,xmax],[modeTempo,modeTempo], color=5
;OPLOT, csi_noMesh, BVfreq_noMesh, linestyle=1, color=5
;OPLOT, csi_noMesh, lamb_freq_noMesh*ldeg*(ldeg+1.0), linestyle=1, color=5



ps_close 




PS_OPEN, "/Users/francesca/Northwestern/Research/create_models_for_tides_and_oscillations/createWDmodelsForOscillationCode_FromMESA_fixedUnitsIn_04022013/output/WD/0p23Msun/log1673/plots/MESAtoCAFein_log1673_0p23Msun_Z0p02_WD_deltaX5eM5in_09222013_thermalOnDynamical.eps",$
          RATIO =0.8, /HICOLOUR
;; ;Load color table
red   = [0,255,195,135, 100, 255, 140, 255,   0,   0,   0,   0,   0,    0, 255, 140, 255, 150]
green = [0,255,195,135, 100, 0  ,   0, 150, 255, 130,   0,   0, 255,  140,   0,   0, 255, 150]
blue  = [0,255,195,135, 100, 0  ,   0,   0,   0,   0, 255, 130, 255,  140, 255, 140,   0,   0]
TVLCT, [red], [green], [blue]
; 0:black, 1:white, 2:light-grey, 3:grey, 4:dark-grey, 5:red, 6:dark- red, 7:orange, 8:green, 9:dark-green
; 10:blue, 11:dark-blue, 12:cyan, 13:dark-cyan, 14:magenta,  15:purple, 16:yellow, 17:dark-yellow
xmin = -0.05
xmax = 1.05
ymin = 1e-4 
ymax = 1e15

PLOT, csi, c4, /ylog, yRange = [ymin,ymax],xRange=[xmin,xmax],YSTYLE=1, XSTYLE=1, Ytitle = '!7s!3!Dth!N/!7s!3!Ddyn!N=c!D4!N',Xtitle = 'r', linestyle = 3

ps_close 





print, "done"

end

