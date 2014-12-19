def uv2intdir(U, V, decl_mag=0, ang_rot=0):
	"""
	#  [int,dir]=uv2intdir(U,V,decl_mag,ang_rot)
	# INVERSO DA
	# FUNCAO INTDIR2UV.M Decompoem o vetor corrente definido pela intensidade
	#                    e direcao (ref. Norte -> Este) considerando a 
	#	             declinacao magnetica e a rotacao do sistema de coordenadas
	# Uso: entre com intensidade, direcao,declinacao magnetica e orientacao do
	# eixo Y, a partir do Norte (0 deg.). Por exemplo, Canal de Sao Sebastiao =
	# 51 deg. A declinacao para oeste e' negativa, p.ex.: -18 deg.  
	# 
	# Se valor de ang_rot nao for suprido, e' assumido que nao ha' rotacao
	# alem disso, se decl_mag tambem nao for suprido nao e' feito correcao
	# magnetica#

	# Author: Roberto Fioravanti Carelli Fontes
	# Depto. de Oceanografia Fisica IOUSP
	Traduzido para python por Rafael Soutelino - IEAPM
	# Laboratorio de Hidrodinamica Costeira (LHICO)

	""" 

	vetor = U + V*1j
	mag = np.abs(vetor)
	dir = np.angle(vetor)
	dir = dir*180/np.pi
	dir = dir - decl_mag + ang_rot

	dir = mod(90-dir, 360)

	return mag, dir


