{ TERMOQUIMICO-R0.PDE

	Este script prevê o campo de temperaturas

   OBS: atenção à compatibilidade das unidades:
   1 cal = 4.1868 J
   1 W = 3.6 kJ/h 

   1) o calor específico - 'ce' - deve ser coerente com a condutibilidade térmica - k - quanto à unidade de energia
   por exemplo, um calor espeífico de 0.2 cal/(g.°C) = 4186.8 * 0.2 = 837.36 J/(kg.°C)
   deve ser compatível com uma condutibilidade térmica de 2.861 W/(m.°C) = 10300 J/(m.h.°C)
   
   2) a condutibilidade térmica - k - deve ser coerente quanto à unidade de tempo
   por exemplo, 10.3 kJ/(m.h.°C) = 2.460 cal/(m.h.°C)

   3) o coeficiente de convecção também deve estar coerente. 
   Por exemplo, um valortípico docoeficiente de convecção é:
   -para concreto-água de cura:	300 W/(m^2.°C)	ou 1080 kJ/(h.m^2.°C)
   -para concreto-ar:						12 W/(m^2.°C)		ou 43.2 kJ/(h.m^2.°C)
   -para concreto-fôrma-ar:			5 W/(m^2.°C)		ou 18.0 kJ/(h.m^2.°C)
	}

TITLE "Calculo Termoquimico"

COORDINATES cartesian2	{ coordinate system, 1D,2D,3D, etc }

VARIABLES	{ system variables }

    Tp												{ Temperatura }
    teq	(threshold=0.002)			{ Tempo equivalente }
    ux
    uy

SELECT			{ method controls }

   terrlim = 1.0e-2
   errlim = 1.0e-4
   painted = on
   CHECK_SS

DEFINITIONS	{ parameter definitions }

   tf = 5*24		{ tempo final de cálculo }

   LX = 1.0
   LY = 1.0

{ DEFINIÇÕES DO PROBLEMA TÉRMICO }

  Text = 25		{ temperatura externa; pode ser especificada uma função senoidal ou outra qualquer do tempo }

   						{ PROPRIEDADES DAS CAMADAS }
   T_lanc=25	{ temperatura de lançamento da camada }
   tz = 0			{ tempo de lançamento }

   						{ PROPRIEDADES TÉRMICAS }
   k = 9.36			{ condutividade térmica }
   rho=2300	{ massa específica }
   ce	=0.9		{ calor específico }
   !Tad				{ temperatura adiabática }
   EaR=4000
   Tad_inf1=16.40021368590
   a1=0.42998953984*24
   c1=1.93329591299
   Tad_inf2=13.59978631410
   a2=1.34507960663*24
   c2=1.57098129317
   Q=rho*ce*(Tad_inf1*teq^c1/(a1^c1+teq^c1) + Tad_inf2*teq^c2/(a2^c2+teq^c2) )

{ DEFINIÇÕES DO PROBLEMA ELÁSTICO }

   E = 20000		{ Young's Modulus - values come later }
   nu	 = 0.2			{ Poisson's Ratio - values come later }
   alpha = 1.e-5	{ Expansion coefficient - values come later }
   Enu = E/(1-nu^2)

{ DEFINIÇÃO DE RELAÇÕES CONSTITUTIVAS PARA MATERIAL ISOTRÓPICO }
{ Expressão do tensor de Hooke para o caso isotrópico }
   C11 = 1
   C22 = 1
   C12 = nu
   C66 = (1-nu)

{ Relações entre deformações e deslocamentos }
   ex = dx(ux)
   ey = dy(uy)
   exy = (dy(ux) + dx(uy))*0.5

{ Lei de Hooke }
   Sxx = Enu*(C11*ex + C12*ey {+ C13*ez} - alpha*(1+nu)*(Tp-T_lanc))
   Syy = Enu*(C12*ex + C22*ey {+ C23*ez} - alpha*(1+nu)*(Tp-T_lanc))
   Txy = Enu*C66*exy

INITIAL VALUES


    Tp = T_lanc 		{ a temperatura inicial será a temperatura de lançamento de cada camada }
    teq= 0.01
    ux  = -1.e-8
    uy  = -1.e-7

EQUATIONS			{ PDE's, one for each variable }

    { Equação do calor }
    Tp: 		div(k*grad(Tp)) + dt(Q) = rho*ce*dt(Tp)

    {Equação da geração interna de calor do processo}
    teq:	dt(teq)=dt(t)*exp(EaR*(1/(298.15)-1/(Tp+273.15)))
then
    { the ux-displacement equation }
    ux:		dx(Sxx) + dy(Txy) = 0.

    { the uy-displacement equation }
    uy:		dx(Txy) + dy(Syy) = 0.

BOUNDARIES       { The domain definition }

   REGION 1       

       start(0,0)	natural(Tp) = 0 value(uy) = 0 line to (LX,0) value(Tp) = Text natural(uy) = 0 value(ux) = 0
							line to (LX,LY)natural(Tp) = 0 nobc(ux)
 							line to (0,LY) value(Tp) = Text value(ux) = 0
							line to close

   TIME from 0.001 TO tf	by 0.01  		{ período de tempo considerado no cálculo }
   !HALT = 0.0001
   !critical tf/2

MONITORS			{ show progress }

PLOTS					{ save result displays }

   for  cycle=5 !t = 0 by tf/100 to tf 

    grid(x+1000*ux,y+1000*uy)
    contour(Tp) as "Campo de temperaturas" !fixed range=(19.0,62.0)
    vector(ux,uy) as "Campo de deslocamentos"
    contour(Sxx) as "Campo de tensao x" !fixed range=(19.0,62.0)
    contour(Syy) as "Campo de tensao y" !fixed range=(19.0,62.0)

    
HISTORIES
    history(Tp) at
            (LX/2,0)
            (LX/4,0)
            (0,0)
    history(ux, uy) at
            (LX/2,0)
            (LX/4,0)
            (0,0)
{    history(Tp, ux, uy) at
            (LX/2,0) export format "#t#r,#i"
    history(Tp, ux, uy) at
            (LX/4,0) export format "#t#r,#i"
    history(Tp, ux, uy) at
            (0,0) export format "#t#r,#i"
  }
END