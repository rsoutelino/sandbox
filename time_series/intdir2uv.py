def intdir2uv(int, dir, decl_mag=0, ang_rot=0):
	"""
	FUNCAO INTDIR2UV Decompoem o vetor corrente definido pela intensidade
	                   e direcao (ref. Norte -> Este) considerando a 
	             declinacao magnetica e a rotacao do sistema de coordenadas
	Uso: entre com intensidade, direcao,declinacao magnetica e orientacao do
	eixo Y, a partir do Norte (0 deg.). Por exemplo, Canal de Sao Sebastiao =
	51 deg. A declinacao para oeste e' negativa, p.ex.: -18 deg. 
	
	Se valor de ang_rot nao for suprido, e' assumido que nao ha' rotacao
	alem disso, se decl_mag tambem nao for suprido nao e' feito correcao
	magnetica
    Author: Roberto Fioravanti Carelli Fontes
    Traduzido para python por Rafael Soutelino - IEAPM
	Depto. de Oceanografia Fisica IOUSP
	Laboratorio de Hidrodinamica Costeira (LHICO)
	"""
	import numpy as np

	dir = dir + decl_mag
	dir = np.mod(dir,360)
	dir = dir*np.pi/180
	
	ang_rot = (ang_rot)*np.pi/180

	u=int*np.sin(dir)
	v=int*np.cos(dir)

	U = u*np.cos(ang_rot) - v*np.sin(ang_rot)
	V = u*np.sin(ang_rot) + v*np.cos(ang_rot)

	return U, V



