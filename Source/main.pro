; 18-09-2013 - 16:57:07
;
; DCE@urLAB
; Copyright (c) 2012,2013, Universidad Politécnica de Madrid
; All rights reserved.
;
; Redistribution and use in source and binary forms, with or without
; modification, are permitted provided that the following conditions are met:
;     * Redistributions of source code must retain the above copyright
;       notice, this list of conditions and the following disclaimer.
;     * Redistributions in binary form must reproduce the above copyright
;       notice, this list of conditions and the following disclaimer in the
;       documentation and/or other materials provided with the distribution.
;     * Neither the name of the copyright holder nor the
;       names of its contributors may be used to endorse or promote products
;       derived from this software without specific prior written permission.
;
; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
; ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
; WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
; DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
; DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
; (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
; ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS

;-------------------------------------------------------------------------------------------

;**************************************************************************************************
;**************************************************************************************************
;<+>

PRO ListOfFunctions

	FORWARD_FUNCTION $
	ANALYZE_AUC,$
	ANALYZE_AUC_ARRAY,$
	ANALYZE_TTM_ARRAY,$
	BRUKERLINE_GETNUMBERVALUE,$
	BRUKERLINE_GETVALUES,$
	CALCULATE_STATISTICS,$
	CHANGE_TRUE,$
	CONVOLUTIONMATRIX,$
	ESTIMATE_CT_ARRAY_TOFTSMODEL,$
	ESTIMATE_CT_ARRAY_YANKEELOVMODEL,$
	ESTIMATE_CT_TOFTSMODEL,$
	ESTIMATE_CT_TOFTSMODEL_CONVOL,$
	ESTIMATE_CT_YANKEELOVMODEL_CONVOL,$
	ESTIMATE_RCE_ARRAY_HOFFMANNMODEL,$
	ESTIMATE_RCE_ARRAY_LARSSONMODEL,$
	ESTIMATE_RCE_HOFFMANNMODEL,$
	ESTIMATE_RCE_LARSSONMODEL,$
	FILTER_GAUSSIAN_1D,$
	FINITE_VALUES,$
	FSC_NORMALIZE,$
	GET_BITUPM_STRARRAY,$
	GET_COLOR,$
	GET_COMPLETE_PATH,$
	GET_DISC,$
	GET_ELAPSEDTIME,$
	GET_FLOATINGFROMSTRING,$
	GET_INTEGERFROMSTRING,$
	GET_LICENSE_ARRAY,$
	GET_NAME,$
	GET_NAME_EXTENSION,$
	GET_NAME_FIELD,$
	GET_OS_INFO,$
	GET_PATH,$
	GET_POINT_YDOWNCONVENTION,$
	GET_POINT_YUPCONVENTION,$
	GET_ROIFROMTXTDATA,$
	GET_WINDOWNUMBER,$
	GETBRUKERHEADER_ARRAY,$
	GETBRUKERHEADER_KEY,$
	GETWIDGET_POSITIONINSCREEN,$
	INTERFACE_LABELINFO_OK,$
	INTERFACE_LABELINFO_CANCEL,$
	GETWIDGET_RELATIVEPOSITION,$
	IDXPOINTS_CHANGEYCONVENTION,$
	IDXPOINTS_FIRSTFROMLOWRES,$
	IDXPOINTSFROMBOX,$
	IDXPOINTSLR_FROM_INDEXES,$
	IMAG2THREECHANNELS,$
	IMAGFROMARRAY,$
	IMIDXS_FROM_INDEXES,$
	INIPOINTFROMINDEXES,$
	INTERFACE_DATAIN_GENERIC,$
	PARAMETRICCP_TO_CURVECP,$
	DCEMRI_READCONFIG,$
	DCEMRI_INITIALPARAMETERS,$
	DCEMRI_CONFIGMODELPARAMETERS,$
	FSTRUCT_CREATE,$
	FWDCEMRI_ACTIVATE_BASES,$
	FWDCEMRI_SUFFIXRESULTNAME,$
	FWDCEMRI_BOXROI_SET,$
	FWDCEMRI_CHECKDRAWROI,$
	FWDCEMRI_DATAFROMSUBIDX,$
	FWDCEMRI_DRAW_IMAGE,$
	FWDCEMRI_DRAW_RCEIMAGE,$
	FWDCEMRI_DYNDATAFROMIDX,$
	FWDCEMRI_EXPORTIMAGE,$
	FWDCEMRI_EXPORTROI,$
	FWDCEMRI_FULLROI_SET,$
	FWDCEMRI_HIDE_IMAGES,$
	FWDCEMRI_LOAD_IMAGE,$
	FWDCEMRI_IMAGEINFO,$
	FWDCEMRI_PLOTSOPENEDTYPE,$
	FWDCEMRI_PARAMETERS_WRITE,$
	FWDCEMRI_PARAMETERS_READ,$
	FWDCEMRI_PARAMETRICAIF_WRITE,$
	FWDCEMRI_PARAMETRICAIF_READ,$
	FWDCEMRI_PLOTROI_SET,$
	FWDCEMRI_PRINT_DATA_3D,$
	FWDCEMRI_PRINT_DATA_2D,$
	FWDCEMRI_REDRAW,$
	FWDCEMRI_DRAWINSAMESIZE,$
	FWDCEMRI_REMAKEINDEXWITHOUTOUTLIERS,$
	FWDCEMRI_ROIKINETICS,$
	FWDCEMRI_ROI_SIZERESOLUTION,$
	FWDCEMRI_RESETRESULTS,$
	FWDCEMRI_RESETROI,$
	FWDCEMRI_RESHAPE_WINDOWS,$
	FWDCEMRI_SENTEVENTS_ADJUSTPALETTE,$
	FWDCEMRI_TEXTWRITE_POINTVAL,$
	FWDCEMRI_DRAWINGEVENTS,$
	FWDCEMRI_PERFORM_MODELANALYSIS,$
	FWBASESTRUCT_WINDOW,$
	FWBASESTRUCT_FRAMES,$
	FWBASESTRUCT_SCALES,$
	FWBASESTRUCT_ZOOM,$
	FWBASESTRUCT_PARAMETERS,$
	INTERFACE_DCEMRI_ABOUT,$
	FWIDGET_PLOT_MODEL_V2,$
	FWDCEMRIRS_ADJUSTSCALES,$
	FWDCEMRIRS_PLOT_MODEL,$
	FWDCEMRIRS_DRAWINGEVENTS_RESULT,$
	FWDCEMRIRS_REDRAW_IM1,$
	FWDCMRI_GETRESULTBASENAME,$
	FWBASESTRUCT_TRANSPARENCY,$
	FWBASESTRUCT_KINETICPARAMETER,$
	INTERFACE_INITIALPARAMETERS,$
	INTERFACE_LABELINFO,$
	GET_BUTTONSICONS,$
	FWPLAYER_DISPLAYIMAGE,$
	FWIDGET_READRAW_EVENT,$
	INTERFACE_READRAW,$
	LOCATE_FIRSTNOVALID_WIDGET,$
	MODEL_CP_BIEXPONENTIAL,$
	MODEL_CP_ORTON2,$
	MODEL_CP_PARKER,$
	MODEL_CP_SCHABEL,$
	MODEL_CP_SCHABELSIMPLIFIED,$
	MODEL_CT_TOFTS,$
	MODEL_CT_TOFTS_CONVOL,$
	MODEL_CT_YANKEELOV_CONVOL,$
	MODEL_RCE_HOFFMANN,$
	MODEL_RCE_LARSSON,$
	MPFIT_CT_TOFTS_BIEXP,$
	MPFIT_CT_TOFTS_CONVOL,$
	MPFIT_CT_TOFTS_HORSFIELD_V1,$
	MPFIT_CT_TOFTS_HORSFIELD_V2,$
	MPFIT_CT_TOFTS_HORSFIELD_V3,$
	MPFIT_CT_YANKEELOV_CONVOL_V1,$
	MPFIT_RCE_HOFFMANN_V1,$
	MPFIT_RCE_HOFFMANN_V2,$
	MPFIT_RCE_HOFFMANN_V3,$
	MPFIT_RCE_HOFFMANN_V4,$
	MPFIT_RCE_HOFFMANN_V5,$
	MPFIT_RCE_LARSSON_V1,$
	NORMALIZEDGAMMA_VARIATE,$
	PLOT_PROFILE,$
	PLOT_PROFILE_DYNMRI,$
	PLOT_PROFILE_KINETICDATA,$
	PLOT_XYOUTS_MATRIX,$
	PLOTPROFILE,$
	PSYM_USERSIM,$
	RCE_SIGNAL_TO_TISSUECONCENTRATION_BARBORIAK,$
	RCE_SIGNALREL_FROM_SIGNAL,$
	READ_AIFFUNCTION,$
	READ_BINARY_BRUKERDATA,$
	READ_BRUKER_T1MAP,$
	READ_BRUKERBIOSPIN_IMAGE,$
	READ_BRUKERBIOSPIN_METHODHEADER,$
	READ_BRUKERBIOSPIN_RECOHEADER,$
	READ_BRUKERBIOSPIN_VISUHEADER,$
	READ_CPMODELCONFIG,$
	READ_DYNAMICDICOM,$
	READ_RAW,$
	READ_RAW_V2,$
	READARRAY_NUMBERS,$
	READFILE_ASCII,$
	READSTRARR_OBTAINVALUES,$
	RELATIVECONTRASTENHANCEMENT,$
	REMOVE_COMMENTS,$
	RESTORE_GRAPHICOBJECTSDATA,$
	ROIMAP_GETPOINT,$
	SETWIDGET_MAXSIZES,$
	SETWIDGET_RELATIVEPOSITION,$
	SIGMOID_CURVE,$
	STANDARDERROR_DIVIDE,$
	STANDARDERROR_MULTIPLY,$
	TRANSLATE_DATA2STRINGS,$
	TRANSLATE_NUMBER2STRING,$
	TRANSLATE_NUMBERARRAY2STRING,$
	UPDATEROIDATA,$
	VIEW_MOSAIC,$
	WRITE_RAW,$
	WRITEFILE_ASCII

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Analyze_AUC, array

	n_points = N_ELEMENTS(array)

	sum_area = 0d

	FOR i=1l, n_points-1 DO BEGIN
		sum_area+=((array[i]+array[i-1])/2) > 0
	ENDFOR

	RETURN, sum_area
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Analyze_AUC_array, data, FRAME_INI=frame_ini, FRAME_END=frame_end, $
			OPT_MAX=opt_max, RELATIVE=relative, ADD_MEAN=add_mean

; analiza el area bajo la curva (desde el comienzo hasta el máximo)

IF KEYWORD_SET(add_mean) THEN opt_add_mean=1l ELSE opt_add_mean=0l
; opción de sumar a los resultados el valor medio en todo el array

CASE SIZE(data, /N_DIMENSIONS) OF
	2 : BEGIN; Array [n_points, n_frames] (100,nframes)
 		type_data = 'ARRAY'
 		n_points = (SIZE(data, /DIMENSIONS))[0]
 		nframes  = (SIZE(data, /DIMENSIONS))[1]
 		im_auc   = FLTARR(n_points+opt_add_mean)
  	END
 	1 : BEGIN
		type_data = 'SIGNAL'
		n_points  = 1
		nframes   = (SIZE(data, /DIMENSIONS))[0]
		im_auc    = FLTARR(1+opt_add_mean)
		data      = REFORM(data, 1, nframes)
	END
 	ELSE : RETURN, -1
ENDCASE

opt_relative = KEYWORD_SET(relative)
IF KEYWORD_SET(opt_max) THEN option_maxvalue=1 ELSE option_maxvalue=0

IF nframes LE frame_ini THEN RETURN, -2

f_ini = frame_ini[0]

IF option_maxvalue EQ 1 THEN BEGIN
	FOR n=0l, n_points-1 DO BEGIN
		arr = REFORM(data[n,f_ini:nframes-1])
		val_max = (WHERE(arr EQ MAX(arr)))[0]
		arr     = arr[0:val_max-1]
		IF opt_relative THEN arr/=arr[0]
		value_auc = Analyze_AUC(arr)
		im_auc[n] = value_auc
	ENDFOR
	IF opt_add_mean THEN BEGIN
		arr_mean = TOTAL(data[*,f_ini:nframes-1],1)/n_points
		val_max  = (WHERE(arr_mean EQ MAX(arr_mean)))[0]
		arr_mean    = arr_mean[0:val_max-1]
		IF opt_relative THEN arr_mean/=arr_mean[0]
		value_auc = Analyze_AUC(arr_mean)
		im_auc[n_points] = value_auc
	ENDIF

ENDIF ELSE BEGIN
	f_end = frame_end[0]
	IF f_end LE f_ini THEN RETURN, -1

	FOR n=0l, n_points-1 DO BEGIN
		arr = REFORM(data[n,f_ini:f_end])
		IF opt_relative THEN arr/=arr[0]
		value_auc = Analyze_AUC(arr)
		im_auc[n] = value_auc
	ENDFOR
	IF opt_add_mean THEN BEGIN
		arr_mean = TOTAL(data[*,f_ini:f_end],1)/n_points
		IF opt_relative THEN arr_mean/=arr_mean[0]
		value_auc = Analyze_AUC(arr_mean)
		im_auc[n_points] = value_auc
	ENDIF

ENDELSE

RETURN , im_auc

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Analyze_TTM_array, data, FRAME_INI=frame_ini, ADD_MEAN=add_mean

; Encuentra el tiempo hasta el frame con valor máximo
; frame_ini : frame de inyección del bolo
; frame_end : Si se


IF KEYWORD_SET(add_mean) THEN opt_add_mean=1l ELSE opt_add_mean=0l
; opción de sumar a los resultados el valor medio en todo el array


CASE SIZE(data, /N_DIMENSIONS) OF
	2 : BEGIN; Array [n_points, n_frames] (100,nframes)
 		type_data = 'ARRAY'
 		n_points = (SIZE(data, /DIMENSIONS))[0]
 		nframes  = (SIZE(data, /DIMENSIONS))[1]
 		im_ttm   = FLTARR(n_points+opt_add_mean)
 	END
 	1 : BEGIN
		type_data = 'SIGNAL'
		n_points  = 1
		nframes   = (SIZE(data, /DIMENSIONS))[0]
		im_ttm    = FLTARR(1+opt_add_mean)
		data      = REFORM(data, 1, nframes)
	END
 	ELSE : RETURN, -1
ENDCASE

opt_relative = KEYWORD_SET(relative)

IF nframes LE frame_ini THEN RETURN, -2

f_ini = frame_ini[0]

FOR n=0l, n_points-1 DO BEGIN
	arr = REFORM(data[n,f_ini:nframes-1])
	Pos_max = (WHERE(arr EQ MAX(arr)))[0]
	im_ttm[n] = pos_max
ENDFOR
IF opt_add_mean THEN BEGIN
	arr_mean = TOTAL(data[*,f_ini:nframes-1],1)/n_points
	Pos_max = (WHERE(arr_mean EQ MAX(arr_mean)))[0]
	im_ttm[n_points] = pos_max
ENDIF

RETURN , im_ttm

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Brukerline_getNumberValue, str_line, NEXT_LINE=next_line, ERROR_ID=error_id

	next_line = 0
	error_id  = 0

	line = STRSPLIT(str_line[0], '=', /EXTRACT)
	IF N_ELEMENTS(line) NE 2 THEN BEGIN
		error_id = -1 & RETURN, -1
	ENDIF
	line = STRTRIM(line[1],2)

	str_pos = STRPOS(line, '(')
	IF str_pos NE -1 THEN BEGIN
		line = STRSPLIT(line, '(', /EXTRACT, /PRESERVE_NULL)
		IF N_ELEMENTS(line) NE 2 THEN  BEGIN
			error_id = -2 & RETURN, -1
		ENDIF
		line = STRTRIM(line[1],2)
		line = STRSPLIT(line, ')', /EXTRACT, /PRESERVE_NULL)
		IF N_ELEMENTS(line) NE 2 THEN  BEGIN
			error_id = -3 & RETURN, -1
		ENDIF
		line = STRTRIM(line[0],2)
		number = FIX(line)
		IF number LE 0 THEN  BEGIN
			error_id = -4 & RETURN, -1
		ENDIF
		next_line = 1
	ENDIF ELSE BEGIN
		number = FLOAT(STRTRIM(line,2))
	    next_line = 0
	ENDELSE


	RETURN, number

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Brukerline_getValues, str_line, n_values, ERROR_ID=error_id

	error_id = 0l

	line = STRSPLIT(str_line[0], ' ', /EXTRACT)
	IF N_ELEMENTS(line) NE n_values THEN BEGIN
		error_id = -1 & RETURN, -1
	ENDIF

	arr_val = FLTARR(n_values)

	FOR i=0l, n_values-1 DO BEGIN
		arr_val[i] = FLOAT(line[i])
	ENDFOR
	IF n_values EQ 1 THEN arr_val = arr_val[0]
	RETURN, arr_val


END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Calculate_Statistics, arr_measurements, arr_model, RANGE=range, STRARR_NAMES=strarr_names

;;http://en.wikipedia.org/wiki/Coefficient_of_determination

; Prueba del estadístico R2, SSE...

; Root mean squared error (RMSE)
; coefficient of determination  (R2 o R-squared)  R2 = 1- SSerr/SStot
; residual sum of squares.  : (SSerr)
; Regression sum of squares : (SSreg)
; Total sum of squares      : (SStot)

; mean error (ME) ; (no absolute values)


IF N_ELEMENTS(range) EQ 0 THEN opt_range = 0 ELSE opt_range=1

n_points = N_ELEMENTS(arr_model)
IF N_ELEMENTS(arr_measurements) NE n_points THEN RETURN, -1


IF opt_range EQ 0 THEN BEGIN

	mean_measurements = MEAN(arr_measurements)

	SSerr = TOTAL((arr_measurements-arr_model)^2)  ;SSerr
	RMSE  = SQRT(SSerr/n_points)

	SStot = TOTAL((arr_measurements  - mean_measurements)^2)
	SSreg = TOTAL((arr_model - mean_measurements)^2)  ; SSreg

	ME = MEAN((arr_measurements-arr_model))

	R_square = 1 - SSerr/SStot

ENDIF ELSE BEGIN
	IF N_ELEMENTS(range) NE 2 THEN RETURN, -1

	n0 = range[0] > 0
	n1 = range[1] < (N_ELEMENTS(arr_model)-1)
	n_points = N_ELEMENTS(arr_model[n0:n1])

	mean_measurements = MEAN(arr_measurements[n0:n1])

	SSerr = TOTAL((arr_measurements[n0:n1]-arr_model[n0:n1])^2)
	RMSE  = SQRT(SSerr/n_points)

	SStot = TOTAL((arr_measurements[n0:n1]  - mean_measurements)^2)
	SSreg = TOTAL((arr_model[n0:n1] - mean_measurements)^2)

	ME = MEAN((arr_measurements[n0:n1]-arr_model[n0:n1]))

	R_square = 1 - SSerr/SStot

ENDELSE

;SStot = SSerr+SSreg
;R_square = SSreg/SStot

strarr_names = ['RMSE', 'ME', 'R_SQUARE', 'SSERR', 'SSREG', 'SSTOT']

RETURN, [RMSE, ME, R_square, SSerr, SSreg, SStot]

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Change_true, imag

; Cambia el parámetro TRUE de una imagen de color
;
; De TRUE=1 a TRUE=3 y viceversa
;
; Útil para manejar imágenes BMP de tres canales


IF N_ELEMENTS(imag) EQ 0 THEN BEGIN & PRINT, 'Datos no válidos' & RETURN, -1 & END

tam = SIZE(imag)

IF (tam[0] NE 3) THEN BEGIN & PRINT, 'Datos no válidos' & RETURN, -1 & END ; No tiene 3 dimensiones
IF (tam[4] NE 1) THEN BEGIN & PRINT, 'Datos no válidos' & RETURN, -1 & END ; No está en formato BYTE

IF (tam[1] EQ 3) THEN BEGIN ; TRUE=1 to 3
	imag_s = BYTARR(tam[2],tam[3],tam[1])
	im_1   = REFORM(imag[0,*,*], tam[2],tam[3])
	im_2   = REFORM(imag[1,*,*], tam[2],tam[3])
	im_3   = REFORM(imag[2,*,*], tam[2],tam[3])
	imag_s = [[[im_1]],[[im_2]],[[im_3]]]
	RETURN, imag_s
ENDIF

IF (tam[3] EQ 3) THEN BEGIN ; TRUE=3 to 1
	imag_s = BYTARR(tam[3],tam[1],tam[2])
	im_1   = REFORM(imag[*,*,0], 1, tam[1],tam[2])
	im_2   = REFORM(imag[*,*,1], 1, tam[1],tam[2])
	im_3   = REFORM(imag[*,*,2], 1, tam[1],tam[2])
	imag_s = [im_1,im_2,im_3]
	RETURN, imag_s
ENDIF

RETURN, -1

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION ConvolutionMatrix, arr_data,  n_cols, SQUARE=square

; Convolution matrix of (n+Np-1) x n (row x columns)
; Np = n_elements(arr_data)

; toeplitz no simétrica

n_points = N_ELEMENTS(arr_data)

n_rows   = n_cols+n_points-1

H = DBLARR(n_cols, n_rows)

; Non symmetric toeptitz matrix (defined by one row and one column)
FOR i=0l, n_cols-1 DO BEGIN
	r_i=i
	r_e=(r_i+n_points-1) < (n_rows-1)
	H[i,r_i:r_e] = arr_data[0:r_e-r_i]
ENDFOR

IF KEYWORD_SET(square) THEN BEGIN
	H=H[*,0:n_cols-1]
ENDIF

RETURN, H

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Estimate_Ct_Array_ToftsModel, data, $
	data_t1map, $
	ST_DATA=st_data, $
	ST_CP_MODEL=st_cp_model,$
	METHOD=method, $
	IMAG_FLAG=imag_flag,  $
	IM_PARAM=im_param, $
	IM_PCERROR=im_pcerror, $
	IM_STATISTICS=im_statistics, $
	ARR_TIME=arr_time, $
	ARR_CP  =arr_cp,   $
	DATA_ERR=data_Err,$
	DATA_CT=data_ct, $
	DATA_ADJUST_CT=data_adjust_ct,$
	ARR_AUX=arr_aux, $
	STAT_RANGE= stat_range,$
	ERROR_ID=error_id, $
	ERROR_STR=error_Str,$
	OPTION_CONVOL=option_convol, $
	COMMON_T10=common_T10,$
	NAMES_STAT  = names_stat,$
	NAMES_PARAM = names_param,$
	ARR_STATUS = arr_status


; Estima los parámetros del modelo de Tofts para el tejido (Ct) en multipixel
; Hereda los métodos de "Estimate_Ct_ToftsModel" pero añade otros para calcular con convolucion "Estimate_Ct_ToftsModel_Convol"

error_id  = 0l
error_str = ''

positive_values = 1

; importante, limitando kep > ktrans
IF N_ELEMENTS(method) EQ 0 THEN method_tofts  = 1 ELSE method_tofts=method[0]

IF KEYWORD_SET(common_t10) EQ 1 THEN opt_common_So = 1 ELSE opt_common_So = 0

opt_convol = KEYWORD_SET(option_convol) ; Con este método se utiliza arr_cp, y no los parametros m1,m2,a1,a2..

n_param     = 5l
n_statistics = 4

ktrans_ini 	= st_data.ktrans_ini
kep_ini    	= st_data.kep_ini
tr         	= st_data.tr
dose       	= st_data.dose   ;  mmol/Kg animal o)
r1         	= st_data.r1     ; (mM s-1)
tau_ini     = st_data.tau_ini
vp_ini      = st_data.vp_ini
angle       = st_data.flip_angle_degrees*!dpi/180.0

; Only in classical Tofts with analytical biexponential Cp
IF opt_convol EQ 0 THEN BEGIN
	IF SIZE(st_cp_model, /TNAME) EQ 'STRUCT' THEN BEGIN
		m1 = st_cp_model.m1
		m2 = st_cp_model.m2
		a1 = st_cp_model.a1
		a2 = st_cp_model.a2
		arr_aux = [a1,a2,m1,m2]
	ENDIF ELSE BEGIN
		error_id = 1 &	error_str = 'Needs biexponential data'
 		RETURN, -1
 	ENDELSE
ENDIF
;------------------------------------------------

; Can work with 2D data (images) or 1D (arrays)

frame_time       = st_data.frame_period            ; minutos
nframe_injection = st_data.nframe_injection      ; donde empieza la inyección del bolo
nframes_adjust   = st_data.nframes_adjust
nframe_end       = (nframe_injection+nframes_adjust-1)


CASE SIZE(data, /N_DIMENSIONS) OF
	2 : BEGIN; Array [n_points, n_frames] (100,nframes)

 		n_points = (SIZE(data, /DIMENSIONS))[0]
 		nframes  = (SIZE(data, /DIMENSIONS))[1]
 		IF SIZE(data_t1map, /N_DIMENSIONS) NE 1 THEN RETURN, -1
 		IF (SIZE(data_t1map, /DIMENSIONS))[0] NE n_points THEN RETURN, -1

		imag_flag   	= BYTARR(n_points+1)
		data_Ct         = FLTARR(n_points+1, nframes)
		Data_adjust_Ct  = FLTARR(n_points+1, nframe_end+1)
		im_pcerror      = FLTARR(n_points+1, n_param)
		im_param 		= FLTARR(n_points+1, n_param) ; ktrans, kep, ve, tau
		im_statistics   = FLTARR(n_points+1, n_statistics)
 	END
 	ELSE : BEGIN
 		error_id = 1 &	error_str = 'Not correct dimensions of data'
 		RETURN, -1
 	END
 ENDCASE

nframe_end     = nframe_end <  (nframes-1)
arr_time       = INDGEN(nframe_end+1)*frame_time
arr_status     = LONARR(n_points)

;-------------------------------------------------------------------------
;; First step, estimates, the mean value in the ROI (to have initial data)

t10_roi    = MEAN(data_t1map)
arr_data   = TOTAL(data,1)/n_points

S_o = MEAN(arr_data[0:nframe_injection])
arr_ct_roi =  RCE_signal_to_tissueConcentration_barboriak(arr_data, T10=t10_roi, TR=tr, R1=r1, S_O=S_o, ANGLE=angle, INJECTION_FRAME=nframe_injection) ; in seconds!!
undefine, arr_data

arr_time_adjust = arr_time[nframe_injection:nframe_end]-(nframe_injection*frame_time)

;--------------------------------------------------------------------------
; Estimation in ROI

CASE opt_convol OF
0 : BEGIN
	result_roi = Estimate_Ct_ToftsModel($
		ARR_VALUES=arr_ct_roi[nframe_injection:nframe_end],$
		ARR_TIME=arr_time_adjust ,$
		KTRANS_INI=ktrans_ini, $
		TAU_INI=tau_ini,$
		PLASMA_VOL_FRACTION=vp_ini,$
		KEP_INI=Kep_ini, $
		VP_INI=vp_ini, $
		DOSE=dose,$
		ARR_AUX=arr_aux,$
		DATA_ERR=data_err,$
		PCERROR=pcerror, $
		QUIET=quiet,$
		STATUS=status,$
		WEIGHTS=arr_weights,$
		METHOD = method_tofts,$
		YFIT=yfit,$
		NAMES_PARAM = names_param_t)
END
1 : BEGIN	; LA AIF se da calculada, no modelada como biexponencial...

	result_roi = Estimate_Ct_ToftsModel_Convol($
		ARR_VALUES=arr_ct_roi[nframe_injection:nframe_end],$
		ARR_TIME = arr_time_adjust,$
		ARR_CP = arr_cp[nframe_injection:nframe_end],  $
		KTRANS_INI=ktrans_ini, $
		KEP_INI=kep_ini, $
		TAU_INI=tau_ini, $
		VP_INI =vp_ini,  $
		QUIET=quiet,  $
		YFIT=yfit,    $
		METHOD=method_tofts, $
		PCERROR=pcerror,  $
		DATA_ERR=data_err,$
		STATUS=status,$
		BESTNORM=bestnorm,$
		NAMES_PARAM = names_param_t)
END
ENDCASE
;-------------------------------------------------------------------------
ktrans_roi = result_roi[0]
kep_roi    = result_roi[1]
ve_roi     = result_roi[2]
tau_roi    = result_roi[3]
vp_roi     = result_roi[4]
;--------------------------------------------------------------------------
data_Ct[n_points,*]        = arr_Ct_roi
Data_adjust_Ct[n_points,nframe_injection:nframe_end] = yfit
im_param[n_points,0:n_param-1]   = result_roi[0:n_param-1]   ; ktrans; kep; ve; tau, vp
im_pcerror[n_points,0:n_param-1] = pcerror[0:n_param-1]
imag_flag[n_points] = status+1
;--------------------------------------------------
im_statistics[n_points,0:1] = (Calculate_Statistics(arr_ct_roi[nframe_injection:nframe_end], yfit, STRARR_NAMES=strarr_names1))[0:1]
im_statistics[n_points,2:3] = (Calculate_Statistics(arr_ct_roi[nframe_injection:nframe_end], yfit, RANGE=stat_range ,STRARR_NAMES=strarr_names2))[0:1]
;--------------------------------------------------

names_stat  = [strarr_names1[0:1], strarr_names2[0:1]+'(Range)']
names_param = names_param_t[0:n_param-1]

;; Segundo paso, hace la estimación punto a punto


FOR n=0l, n_points-1 DO BEGIN
	arr_data = REFORM(data[n, *])
	t10_point = data_t1map[n]
	;----------------------------
	IF t10_point NE 0 THEN BEGIN
		IF opt_common_So EQ 0 THEN S_o = MEAN(arr_data[0:nframe_injection]) ; is not calculated when is common to all points

		arr_Ct = RCE_signal_to_tissueConcentration_barboriak(arr_data, T10=t10_point, TR=tr, R1 = r1, S_o=S_o, ANGLE=angle, INJECTION_FRAME=nframe_injection) ; Hacer esto ahora con opción de flip angle

		CASE opt_convol OF
		0 : BEGIN
			;-----------------------------------------------------------------------------------------------
			result = Estimate_Ct_ToftsModel(ARR_VALUES=arr_ct[nframe_injection:nframe_end],$
				ARR_TIME=arr_time_adjust,PLASMA_VOL_FRACTION=vp_ini, $
				KTRANS_INI=Ktrans_roi, KEP_INI=Kep_roi, VP_INI=vp_roi, TAU_INI=tau_roi,$
				DOSE=dose, ARR_AUX=arr_aux, YFIT=yfit, DATA_ERR=data_err,$
				QUIET=1, PCERROR=pcerror, STATUS=status, METHOD=method_tofts)
			;-------------------------------------------------
		END
		1 : BEGIN
			result = Estimate_Ct_ToftsModel_convol(ARR_VALUES=arr_ct[nframe_injection:nframe_end],$
				ARR_CP = arr_cp[nframe_injection:nframe_end], ARR_TIME=arr_time_adjust, $
				KTRANS_INI=Ktrans_roi, KEP_INI=Kep_roi, VP_INI=vp_roi, TAU_INI=tau_roi, $
				YFIT=yfit, DATA_ERR=data_err,$
				QUIET=1, PCERROR=pcerror, STATUS=status, METHOD=method_tofts)
			;-------------------------------------------------
		END
		ENDCASE
		;arr_ct_estimated = TEMPORARY(yfit)
		;arr_ct_estimated = Model_Ct_Tofts(DOSE=dose, KTRANS=result[0], KEP=result[1],$
		;	ARR_AUX=arr_aux, ARR_TIME=arr_time_adjust)
		;-------------------------------------------------
		data_Ct[n,*]        = arr_Ct
		Data_adjust_Ct[n,nframe_injection:nframe_end] = yfit
		IF N_ELEMENTS(result) NE 1 THEN BEGIN
			im_param[n,0:n_param-1]     = result[0:n_param-1]   ; ktrans; kep; ve; tau, vp
		ENDIF
		im_pcerror[n,0:n_param-1] = pcerror[0:n_param-1]
		imag_flag[n] = status+1
		;--------------------------------------------------
		im_statistics[n,0:1] = (Calculate_Statistics(arr_ct[nframe_injection:nframe_end], yfit, STRARR_NAMES=strarr_names1))[0:1]
		im_statistics[n,2:3] = (Calculate_StatiStics(arr_ct[nframe_injection:nframe_end], yfit, RANGE=stat_range,STRARR_NAMES=strarr_names2))[0:1]
		;--------------------------------------------------
		arr_status[n] = TEMPORARY(status)
	ENDIF
ENDFOR


RETURN, 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>


FUNCTION Estimate_Ct_Array_YankeelovModel, data, $
	data_t1map, $
	ST_DATA=st_data, $
	ST_CP_MODEL=st_cp_model,$
	IMAG_FLAG=imag_flag,  $
	IM_PARAM=im_param, $
	ARR_TIME=arr_time, $
	IM_PCERROR=im_pcerror, $
	DATA_ERR=data_Err,$
	DATA_CT_T=data_ct_t, $
	DATA_ADJUST_CT_T=data_adjust_ct_t,$
	ARR_CT_RR=arr_ct_rr,$
	ARR_AUX=arr_aux, $
	STAT_RANGE= stat_range,$
	IM_STATISTICS=im_statistics, $
	ERROR_ID=error_id, $
	ERROR_STR=error_Str,$
	OPTION_CONVOL=option_convol, $
	COMMON_T10=common_T10,$
	ARR_SIGNAL_RR = arr_signal_rr, $ ; la ROI de la RR se ha procesado antes (quitamos complicacion aqui
	KTRANS_RR = ktrans_rr,$
	KEP_RR    = kep_rr,$
	NAMES_PARAM = names_param,$
	NAMES_STAT  = names_stat


; Estima los parámetros del modelo de Yankeelov (Reference model) para el tejido (Ct) en multipixel
; Hereda los métodos de "Estimate_Ct_YankeelovModel_convol"

error_id  = 0l
error_str = ''
positive_values = 1

IF KEYWORD_SET(common_t10) EQ 1 THEN opt_common_So = 1 ELSE opt_common_So = 0

n_param     = 5l
n_statistics = 4

ktrans_ini 	= st_data.ktrans_ini
kep_ini    	= st_data.kep_ini
tr         	= st_data.tr
dose       	= st_data.dose   ;  mmol/Kg animal o)
r1         	= st_data.r1     ; (mM s-1)
tau_ini     = st_data.tau_ini
vp_ini      = st_data.vp_ini
angle       = st_data.flip_angle_degrees*!dpi/180.0
;------------------------------------------------
; Can work with 2D data (images) or 1D (arrays)

frame_time       = st_data.frame_period            ; minutos
nframe_injection = st_data.nframe_injection      ; donde empieza la inyección del bolo
nframes_adjust   = st_data.nframes_adjust
nframe_end       = (nframe_injection+nframes_adjust)

CASE SIZE(data, /N_DIMENSIONS) OF
	2 : BEGIN; Array [n_points, n_frames] (100,nframes)

 		n_points = (SIZE(data, /DIMENSIONS))[0]
 		nframes  = (SIZE(data, /DIMENSIONS))[1]
 		IF SIZE(data_t1map, /N_DIMENSIONS) NE 1 THEN RETURN, -1
 		IF (SIZE(data_t1map, /DIMENSIONS))[0] NE n_points THEN RETURN, -1

		imag_flag   	= BYTARR(n_points+1)
		data_Ct_t       = FLTARR(n_points+1, nframes)
		Data_adjust_Ct  = FLTARR(n_points+1, nframe_end+1)
		im_pcerror      = FLTARR(n_points+1, n_param)
		im_param 		= FLTARR(n_points+1, n_param) ; ktrans, kep, ve, tau
		im_statistics      = FLTARR(n_points+1, n_statistics)
 	END
 	ELSE : BEGIN
 		error_id = 1 &	error_str = 'Not correct dimensions of data'
 		RETURN, -1
 	END
 ENDCASE

nframe_end     = nframe_end <  (nframes-1)
arr_time       = INDGEN(nframe_end+1)*frame_time

;-------------------------------------------------------------------------
;; First step, estimates, the mean value in the ROI (to have initial data)

t10_roi    = MEAN(data_t1map)
arr_data   = TOTAL(data,1)/n_points

S_o = MEAN(arr_data[0:nframe_injection])
arr_ct_roi =  RCE_signal_to_tissueConcentration_barboriak(arr_data, T10=t10_roi, TR=tr, R1=r1, S_O=S_o, ANGLE=angle, INJECTION_FRAME=nframe_injection) ; in seconds!!
undefine, arr_data

; also translate signal in reference region to Ct_RR
S_o_rr     = MEAN(arr_signal_rr[0:nframe_injection])
arr_ct_rr  =  RCE_signal_to_tissueConcentration_barboriak(arr_signal_rr, T10=t10_roi, TR=tr, R1=r1, S_O=S_o_rr, ANGLE=angle, INJECTION_FRAME=nframe_injection) ; in seconds!!

arr_time_adjust = arr_time[nframe_injection:nframe_end]-(nframe_injection*frame_time)

;--------------------------------------------------------------------------
; Estimation in ROI

result_roi = Estimate_Ct_YankeelovModel_Convol($
		ARR_VALUES=arr_ct_roi[nframe_injection:nframe_end],$
		ARR_TIME = arr_time_adjust,$
		ARR_CT_RR  = arr_ct_rr[nframe_injection:nframe_end],  $
		KTRANS_RR = ktrans_rr,$
		KEP_RR = kep_rr, $
		KTRANS_INI=ktrans_ini, $
		KEP_INI=kep_ini, $
		TAU_INI=tau_ini, $
		VP_INI =vp_ini,  $
		QUIET=quiet,  $
		YFIT=yfit,    $
		METHOD=method_tofts, $
		PCERROR=pcerror,  $
		DATA_ERR=data_err,$
		STATUS=status,$
		BESTNORM=bestnorm,$
		NAMES_PARAM = names_param_t)

;-------------------------------------------------------------------------
ktrans_roi = result_roi[0]
kep_roi    = result_roi[1]
ve_roi     = result_roi[2]
tau_roi    = result_roi[3]
vp_roi     = result_roi[4]
;--------------------------------------------------------------------------
data_Ct_t[n_points,*]       = arr_Ct_roi
Data_adjust_Ct[n_points,nframe_injection:nframe_end] = yfit
im_param[n_points,0:n_param-1]     = result_roi[0:n_param-1]   ; ktrans; kep; ve; tau, vp
im_pcerror[n_points,0:n_param-1]   = pcerror[0:n_param-1]
imag_flag[n_points] = status+1
;--------------------------------------------------
im_statistics[n_points,0:1] = (Calculate_Statistics(arr_ct_roi[nframe_injection:nframe_end], yfit, STRARR_NAMES=strarr_names1))[0:1]
im_statistics[n_points,2:3] = (Calculate_Statistics(arr_ct_roi[nframe_injection:nframe_end], yfit, RANGE=stat_range, STRARR_NAMES=strarr_names2))[0:1]
;--------------------------------------------------

names_stat  = [strarr_names1[0:1], strarr_names2[0:1]+'(Range)']
names_param = names_param_t[0:n_param-1]

;; Segundo paso, hace la estimación punto a punto


FOR n=0l, n_points-1 DO BEGIN
	arr_data = REFORM(data[n, *])
	t10_point = data_t1map[n]
	;----------------------------
	IF t10_point NE 0 THEN BEGIN
		IF opt_common_So EQ 0 THEN S_o = MEAN(arr_data[0:nframe_injection]) ; is not calculated when is common to all points

		arr_Ct = RCE_signal_to_tissueConcentration_barboriak(arr_data, T10=t10_point, TR=tr, R1 = r1, S_o=S_o, ANGLE=angle, INJECTION_FRAME=nframe_injection)
		 ; Hacer esto ahora con opción de flip angle

		result = Estimate_Ct_YankeelovModel_convol(ARR_VALUES=arr_ct[nframe_injection:nframe_end],$
				ARR_CT_RR = arr_ct_RR[nframe_injection:nframe_end], ARR_TIME=arr_time_adjust, $
				KTRANS_INI=Ktrans_roi, KEP_INI=Kep_roi, VP_INI=vp_roi, TAU_INI=tau_roi, $
				KTRANS_RR = ktrans_rr, KEP_RR = kep_rr, YFIT=yfit, DATA_ERR=data_err,$
				QUIET=1, PCERROR=pcerror, STATUS=status, METHOD=method_tofts)
		;----------------------------------------------------------------------------------

		;arr_ct_estimated = TEMPORARY(yfit)
		;arr_ct_estimated = Model_Ct_Tofts(DOSE=dose, KTRANS=result[0], KEP=result[1],$
		;	ARR_AUX=arr_aux, ARR_TIME=arr_time_adjust)
		;----------------------------------------------------------------------------------
		data_Ct_t[n,*]        = arr_Ct
		Data_adjust_Ct[n,nframe_injection:nframe_end] = yfit
		IF N_ELEMENTS(result) NE 1 THEN BEGIN
			im_param[n,0:n_param-1] = result[0:n_param-1]   ; ktrans; kep; ve; tau, vp
		ENDIF
		im_pcerror[n,0:n_param-1] = pcerror[0:n_param-1]
		imag_flag[n] = status+1
		;----------------------------------------------------------------------------------
		im_statistics[n,0:1] = (Calculate_Statistics(arr_ct[nframe_injection:nframe_end], yfit))[0:1]
		im_statistics[n,2:3] = (Calculate_Statistics(arr_ct[nframe_injection:nframe_end], yfit, RANGE=stat_range))[0:1]
		;----------------------------------------------------------------------------------
	ENDIF
ENDFOR

RETURN, 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Estimate_Ct_ToftsModel,$
		ARR_VALUES=arr_values,  $
		ARR_TIME=arr_time, $
		DOSE=dose, $
		ARR_AUX=arr_aux, $
		ARR_INI=arr_ini, $
		QUIET=quiet,  $
		DATA_ERR=data_err, $
		KTRANS_INI=ktrans_ini, $
		KEP_INI=kep_ini, $
		VP_INI =vp_ini,  $
		TAU_INI=tau_ini, $
		PCERROR=pcerror, $
		BESTNORM=bestnorm,  $
		STATUS=status,  $
		WEIGHTS=weights,  $
		METHOD=method,$
		PLASMA_VOL_FRACTION=plasma_vol_fraction, $
		YFIT=yfit,$
		NAMES_PARAM=names_param


; Estima los valores de Ktrans con el modelo de tofts,para concentracion en el tejido (c_t)

; Si PLASMA_VOL_FRACTION se proporciona, estima el modelo modificado de Tofts con
;        el término adicional Vp*Cp, pero Vp conocido

; Con method=1
;     Estimación estándar de Tofts

; Con method=2
;     Estimación de Tofts modificada (término adicional Vp*Cp)) Vp desconocido

; Con Method 3
; - Parámetro adicional (tau shift) de retardo de Cp biexponencial (Orton et al, 2007)
;
; Con method 4
;  -Estimación de Tofts modificada (término adicional Vp*Cp))  Vp desconocido, + tau_shift
;
; -PLASMA_VOL_FRACTION = tanto por uno (entre 0 y 1) del volumen de plasma. Necesario para el método 5
;

; -ARR_AUX - [[a1,a2,m1,m2] ; Parámetros del modelo biexponencial de la AIF
; -DOSE    - Dosis inyectada, en mmol/kg del animal


IF N_ELEMENTS(method) EQ 0 THEN opt_method=1 ELSE opt_method=FIX(method[0])


IF N_ELEMENTS(arr_values) NE N_ELEMENTS(arr_time) THEN RETURN, -1


IF MAX(opt_method EQ [1,2,3,4]) NE 1 THEN BEGIN
	PRINT, 'Not valid method' &	RETURN, -1
ENDIF

n_frames = N_ELEMENTS(arr_values)

IF N_ELEMENTS(data_err) NE 0 THEN $
	rerr = REPLICATE(data_err[0], n_frames)


ktrans_ini = ktrans_ini[0]
kep_ini    = kep_ini[0]
vp_ini     = vp_ini[0]
tau_ini    = tau_ini[0]
ve_ini     = ktrans_ini/kep_ini;

IF N_ELEMENTS(plasma_vol_fraction) EQ  0 THEN vp = 0.0 ELSE vp =  plasma_vol_fraction[0]

IF kep_ini LE ktrans_ini THEN RETURN, -1

X = [Dose[0], arr_aux[0:3], arr_time]
R = arr_values


result_total = DBLARR(5) ; 5  valores:
perror_total = DBLARR(5) ; 5  valores de l-sigma error: [ktrans, kep, ve, tau, vp] ;(tau=0, vp=0, cuando no se calcula)

names_param = ['KTRANS', 'KEP', 'VE', 'VP', 'TAU']
;--------------------------------------------------------------------
ve_ini = ktrans_ini/kep_ini;
Start   =DOUBLE([ktrans_ini, ve_ini, vp_ini, tau_ini])

parinfo =REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0]}, 4)
parinfo[*].value = start
; los dos parametros limitados cantidades positivas

parinfo[1].limited = [1,1]; ve limited to [0,1]
parinfo[1].limits = [0,1D];

parinfo[2].limited = [1,1]; vp limited to [0,1]
parinfo[2].limits = [0,1D];

;-------------------------------------------------------------------------
CASE method OF
1 : fixed_parameters = 2
2 : fixed_parameters = 1
3 : fixed_parameters = 1
4 : fixed_parameters = 0
END
;-------------------------------------------------------------------------

IF (method EQ 1) OR (method EQ 3) THEN BEGIN
		parinfo[2].fixed = 1 ;
		parinfo[2].value = vp;  Tofts stándar (vp=0)
ENDIF
IF (method EQ 1) OR (method EQ 2) THEN BEGIN
		parinfo[3].fixed = 1;
		parinfo[3].value = 0;  Tau = 0 no se estima, Tofts stándar (vp=0)
ENDIF
;-------------------------------------------------------------------------


result = MPFITFUN('MPFIT_CT_TOFTS_BIEXP', X, R, rerr, start, QUIET=quiet,$
	COVAR=covar, STATUS=status, FTOL=1d-20, MPSIDE=mpside, WEIGHTS=weights,PARINFO=parinfo,$
	BESTNORM=bestnorm, DOF=dof, PERROR=perror, YFIT=yfit)

IF N_ELEMENTS(result) EQ 4 THEN BEGIN
	result_total[0] = result[0] ;Ktrans_estimated, kep_estimated
	result_total[1] = result[0]/result[1] ;kep     = Ktrans_estimated/vep_estimated
	result_total[2] = result[1]
	result_total[3:4]  =result[2:3]
ENDIF
IF N_ELEMENTS(perror) EQ 4 THEN BEGIN
	perror_total[0] = perror[0]
	perror_total[1] = StandardError_divide(MEANS=DOUBLE(result[0:1]), SD=DOUBLE(perror[0:1]))
	perror_total[2] = perror[1]
	perror_total[3:4] = perror[2:3]
ENDIF



degrees_of_freedom = N_ELEMENTS(R)-N_ELEMENTS(start)
pcerror = perror_total*SQRT(bestnorm/degrees_of_freedom)

RETURN, result_total ; ktrans, kep, ve, (y opcionalmente tau) vp
;5  valores: [ktrans, kep, ve, tau, vp]

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Estimate_Ct_ToftsModel_Convol,$
	ARR_VALUES=arr_values,$
	ARR_TIME = arr_time,$
	ARR_CP = arr_cp,  $
	KTRANS_INI=ktrans_ini, $
	KEP_INI=kep_ini,  $
	VP_INI =vp_ini,   $
	TAU_INI=tau_ini,  $
	PLASMA_VOL_FRACTION=plasma_vol_fraction, $
	QUIET=quiet,      $
	YFIT=yfit,        $
	METHOD=method,    $
	INTERPOLATION=interpolation,$
	PCERROR=pcerror,  $
	DATA_ERR=data_err,$
	STATUS=status,    $
	BESTNORM=bestnorm,$
	NAMES_PARAM=names_param

; Estimación del modelo de Tofts con AIF genérica, calculando la convolución
;
; Con Method 1
;   Estimación estándar de Tofts
;
; Con method=2
;   Estimación de Tofts modificada (término adicional Vp*Cp)) Vp desconocido

; Con Method 3
; - Parámetro adicional (tau shift) de retardo de Cp biexponencial (Orton et al, 2007)
;
; Con method 4
;  -Estimación de Tofts modificada (término adicional Vp*Cp)) Vp desconocido, + tau_shift
;
; Con method 5
;  -Método de Horsfield y Morgan (2004)


IF N_ELEMENTS(method) EQ 0 THEN opt_method=1 ELSE opt_method= FIX(method[0])

valid_methods = [1,2,3,4,5]
IF MAX(method EQ valid_methods) NE 1 THEN RETURN,-1


IF N_ELEMENTS(arr_values) NE N_ELEMENTS(arr_time) THEN RETURN, -1
IF N_ELEMENTS(arr_values) NE N_ELEMENTS(arr_cp) THEN RETURN, -1


ktrans_ini = ktrans_ini[0]
kep_ini    = kep_ini[0]
vp_ini     = vp_ini[0]
tau_ini    = tau_ini[0]
ve_ini     = ktrans_ini/kep_ini;

result_total = DBLARR(5) ; 5  valores: [ktrans, kep, ve, tau, vp] (tau no se calcula pero se pone por compatibilidad.
perror_total = DBLARR(5)

IF N_ELEMENTS(plasma_vol_fraction) EQ  0 THEN vp = 0.0 ELSE vp =  plasma_vol_fraction[0]

names_param = ['KTRANS', 'KEP', 'VE', 'VP', 'TAU']
;-------------------------------------------------------------------------
CASE opt_method OF
1 : fixed_parameters = 2
2 : fixed_parameters = 1
3 : fixed_parameters = 1
4 : fixed_parameters = 0
5 : fixed_parameters = 2
END
;-------------------------------------------------------------------------

IF opt_method LE 4 THEN BEGIN

	X = [REFORM(arr_time), REFORM(arr_cp)]
	R = arr_values ; arr_sampled_Ct

	ve_ini = ktrans_ini/kep_ini;
	Start   =DOUBLE([ktrans_ini, ve_ini, vp_ini, tau_ini])

	parinfo =REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0]}, 4)
	parinfo[*].value = start
	; los dos parametros limitados cantidades positivas

	parinfo[1].limited = [1,1]; ve limited to [0,1]
	parinfo[1].limits = [0,1D];

	parinfo[2].limited = [1,1]; vp limited to [0,1]
	parinfo[2].limits = [0,1D];


	IF (opt_method EQ 1) OR (opt_method EQ 3) THEN BEGIN
			parinfo[2].fixed = 1 ;
			parinfo[2].value = vp;  Tofts stándar (vp=0)
	ENDIF
	IF (opt_method EQ 1) OR (opt_method EQ 2) THEN BEGIN
			parinfo[3].fixed = 1;
			parinfo[3].value = 0;  Tau = 0 no se estima, Tofts stándar (vp=0)
	ENDIF
	;-------------------------------------------------------------------------

	result = MPFITFUN('MPFIT_CT_TOFTS_CONVOL', X, R, rerr, start, QUIET=quiet, YFIT=yfit,$
		BESTNORM=bestnorm, DOF=dof, PERROR=perror, PARINFO=parinfo, STATUS=status)

	arr_Ct_est  = Model_Ct_Tofts_convol(ARR_TIME=arr_time, ARR_CP=arr_cp, KTRANS=result[0], KEP=result[0]/result[1], NORMALIZE=0, METHOD=1)

	;IF yfit[N_ELEMENTS(yfit)-1] EQ 0 THEN BEGIN & ok = DIALOG_MESSAGE('error temporal') & RETURN,-1 & ENDIF

	IF N_ELEMENTS(result) EQ 4 THEN BEGIN
		result_total[0] =result[0]           ;Ktrans_estimated, kep_estimated
		result_total[1] =result[0]/result[1] ;kep_estimated, Ktrans_estimated/ve_estimated
		result_total[2:4] =result[1:3]       ;ve_estimated, vp and tau
	ENDIF
	IF N_ELEMENTS(perror) EQ 4 THEN BEGIN
		perror_total[0] = perror[0]
		perror_total[1] = StandardError_divide(MEANS=DOUBLE(result[0:1]), SD=DOUBLE(perror[0:1]))
		perror_total[2:4] = perror[1:3]
	ENDIF

	degrees_of_freedom = N_ELEMENTS(R)-N_ELEMENTS(start)-fixed_parameters
	pcerror = perror_total*SQRT(bestnorm/degrees_of_freedom)

ENDIF
;-------------------------------------------------------------------------
IF opt_method EQ 5 THEN BEGIN

	X = [REFORM(arr_values), REFORM(arr_cp)]
	R = arr_values ; arr_sampled_Ct

	inc_time = (arr_time[N_ELEMENTS(arr_time)-1]-arr_time[0])/(N_ELEMENTS(arr_time)-1)

	ktrans_dt_ini = ktrans_ini*inc_time

	Eparam_ini    = EXP(-ktrans_dt_ini/ve_ini)
	Start   = DOUBLE([ktrans_dt_ini, ve_ini])

	parinfo =REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0]}, 2)
	parinfo[*].value = start
	; los dos parametros limitados cantidades positivas
	fixed_parameters = 0

	parinfo[1].limited = [1,1]; exp(negative_number)
	parinfo[1].limits = [0,1D];

	IF N_ELEMENTS(interpolation) EQ 0 THEN opt_interpol=1 ELSE opt_interpol = interpolation[0]
	CASE opt_interpol OF
	1 : function_mpfit = 'MPFIT_CT_TOFTS_HORSFIELD_V1'
	2 : function_mpfit = 'MPFIT_CT_TOFTS_HORSFIELD_V2'
	3 : function_mpfit = 'MPFIT_CT_TOFTS_HORSFIELD_V3'
	ENDCASE

	result = MPFITFUN( function_mpfit, X, R, rerr, start, QUIET=quiet, YFIT=yfit, MPSIDE=2, $
		BESTNORM=bestnorm, DOF=dof, PERROR=perror, PARINFO=parinfo, STATUS=status)

	IF N_ELEMENTS(result) EQ 2 THEN BEGIN
		ktrans_dt = result[0]
		ve        = result[1]

		ktrans =  ktrans_dt/inc_time
		kep    =  ktrans/ve

		result_total[0:2] = [ktrans, kep, ve]

		ktrans_stdev = perror[0]/inc_time ;standard error divided by a number
		ve_stdev = perror[1]

		perror_total[0] = ktrans_stdev
		perror_total[1] = StandardError_divide(MEANS=DOUBLE([ktrans, ve]), SD=DOUBLE([ktrans_stdev,ve_stdev]))
		perror_total[1] = ve_stdev
		; ????

	ENDIF

	degrees_of_freedom = N_ELEMENTS(R)-N_ELEMENTS(start)-fixed_parameters
	pcerror = perror_total*SQRT(bestnorm/degrees_of_freedom)


ENDIF


RETURN, result_total


END


;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Estimate_Ct_YankeelovModel_Convol,$
	ARR_VALUES = arr_values,  $ ; like arr_ct
	ARR_TIME = arr_time,  $
	ARR_CT_RR= arr_ct_rr, $
	KTRANS_RR= ktrans_rr, $
	KEP_RR   = kep_rr,    $
	KTRANS_INI=ktrans_ini,$
	KEP_INI=kep_ini,  $
	VP_INI =vp_ini,   $
	TAU_INI=tau_ini,  $
	QUIET=quiet,      $
	YFIT=yfit,        $
	METHOD=method,    $
	PCERROR=pcerror,  $
	DATA_ERR=data_err,$
	STATUS=status,    $
	BESTNORM=bestnorm,$
	NAMES_PARAM=names_param

; Estimación del modelo de RR de yankeelov (arr_ct_rr genérica), calculando la convolución

IF N_ELEMENTS(method) EQ 0 THEN opt_method=1 ELSE opt_method= FIX(method[0])


IF N_ELEMENTS(arr_values) NE N_ELEMENTS(arr_time)  THEN RETURN, -1
IF N_ELEMENTS(arr_ct_rr) NE N_ELEMENTS(arr_time) THEN RETURN, -1

ve_rr = ktrans_rr[0]/kep_rr[0]

X = [ktrans_rr[0], ve_rr[0], REFORM(arr_time), REFORM(arr_ct_rr)]
R = arr_values ; arr_sampled_Cp

ktrans_ini = ktrans_ini[0]
kep_ini    = kep_ini[0]
vp_ini     = vp_ini[0]
tau_ini    = tau_ini[0]
ve_ini     = ktrans_ini/kep_ini;

result_total = DBLARR(5) ; 5  valores: [ktrans, kep, ve, tau, vp] (tau no se calcula pero se pone por compatibilidad.
perror_total = DBLARR(5)

names_param = ['KTRANS', 'KEP', 'VE', 'TAU', 'VP']

CASE opt_method OF
1 : BEGIN
	Start   =DOUBLE([ktrans_ini, ve_ini])
	parinfo =REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0]}, 2)
	parinfo[*].value = start
	; los dos parametros limitados cantidades positivas
	fixed_parameters = 0

	result = MPFITFUN('MPFIT_CT_YANKEELOV_CONVOL_V1', X, R, rerr, start, QUIET=quiet, YFIT=yfit,$
		BESTNORM=bestnorm, DOF=dof, PERROR=perror, PARINFO=parinfo, STATUS=status)

	IF N_ELEMENTS(result) EQ 2 THEN BEGIN
		result_total[0] =result[0]           ;Ktrans_estimated, kep_estimated
		result_total[1] =result[0]/result[1] ;kep_estimated, Ktrans_estimated/ve_estimated
		result_total[2] =result[1]           ;ve_estimated
		result_total[3] =0  ; tau
		result_total[4] =0  ; vp
	ENDIF
	IF N_ELEMENTS(perror) EQ 2 THEN BEGIN
		perror_total[0] = perror[0]
		perror_total[1] = StandardError_divide(MEANS=DOUBLE(result[0:1]), SD=DOUBLE(perror[0:1]))
		perror_total[2] = perror[1]
		perror_total[3] = 0
		perror_total[4] = 0
	ENDIF

END
ENDCASE

degrees_of_freedom = N_ELEMENTS(R)-N_ELEMENTS(start)-fixed_parameters
pcerror = perror_total*SQRT(bestnorm/degrees_of_freedom)


RETURN, result_total


END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Estimate_RCE_Array_HoffmannModel, data, $
	ARR_TIME=arr_time,    $
	ARR_AIF =arr_aif,     $
	ST_DATA=st_data,      $
	DATA_ERR=data_err,    $
	IM_PARAM=im_param,     $
	IM_PCERROR=im_pcerror, $

	STAT_RANGE= stat_range,$
	IM_STATISTICS= im_statistics,$

	NO_ROI_ESTIMATION=no_roi_estimation,$
	DATA_SIGNAL_RCE=data_signal_rce, $
	DATA_ADJUST=data_adjust, $
	METHOD=method,$
	QUIET=quiet,  $
	TAU_INFUSION=tau_infusion,$
	ERROR_ID=error_id, ERROR_STR=error_str,$
	NAMES_PARAM=names_param,$
	NAMES_STAT =names_stat




positive_values = 1

st_data_hoffmann = {$
	kep_ini :st_data.keph_ini,$
	kel_ini :st_data.kelh_ini,$
	ah_ini  :st_data.Ah_ini, $
	bh_ini  :st_data.Bh_ini, $
	tau_ini :st_data.tau_ini,$
	tau_infusion:st_data.tau_infusion}


frame_time       = st_data.frame_period          ; minutos
nframe_injection = st_data.nframe_injection      ; donde empieza la inyección del bolo
nframes_adjust   = st_data.nframes_adjust
nframe_end       = (nframe_injection+nframes_adjust)



; data_signal_RCE son valores de contraste relativo en tanto por uno (ya deben de estar normalizados)
; El primer valor no tiene porqué ser uno.. (está en media con respecto a antes del contraste)

error_id = 0l & error_str = ''

IF KEYWORD_SET(method) THEN opt_method=FIX(method[0]) ELSE  opt_method = 1l

; Data_err es un valor estimado de DESVIACIÓN ESTÁNDAR de los datos..( para entrada ERR de mpfitfun)
; estima el modelo de Hoffmann para varios puntos
; data_signal_RCE en % por uno

n_parameters = 6l
n_statistics = 4l ; RMSE,ME,  RMSE(range),ME (range)

CASE SIZE(REFORM(data), /N_DIMENSIONS) OF
	2 : BEGIN; Array [n_points, n_frames] (100,nframes)
 		type_data       = 'ARRAY'
 		n_points        = (SIZE(REFORM(data), /DIMENSIONS))[0]
 		nframes         = (SIZE(REFORM(data), /DIMENSIONS))[1]
 		im_param        = FLTARR(n_points+1,n_parameters)   ; +1 dimensions to include total ROI
 		im_pcerror      = FLTARR(n_points+1,n_parameters)
 		im_statistics   = FLTARR(n_points+1,n_statistics)
 		data_signal_RCE = FLTARR(N_points+1,nframes)
		data_adjust     = FLTARR(n_points+1,nframes)
		data_adjust[*,0:nframe_injection]=1.0
	END
 	ELSE : BEGIN
 		error_id = 1 &	error_str = 'Not correct dimensions of data'
 		RETURN, -1
 	END
ENDCASE
nframe_end      = nframe_end  <  (nframes-1)
arr_time        = INDGEN(nframes)*1.0*frame_time
arr_time_adjust = arr_time[nframe_injection:nframe_end]-(frame_time*nframe_injection)


Data_signal_RCE[0:n_points-1,*]  = RCE_signalRel_from_Signal(data, FRAME_INJECTION=nframe_injection) ; only calculates RCE in % per one
Data_signal_RCE[n_points,*]      = RCE_signalRel_from_Signal(TOTAL(data,1)/n_points, FRAME_INJECTION=nframe_injection) ; mean value

; data_signal_RCE son valores de contraste relativo en tanto por uno (ya deben de estar normalizados)
; El primer valor no tiene porqué ser uno.. (está en media con respecto a antes del contraste)

;------------------------------------------------------------------------------------
IF KEYWORD_SET(no_roi_estimation) EQ 0 THEN BEGIN
	; Primer paso, estima los parámetros medios para todo el ROI

	arr_RCE_mean = REFORM(Data_signal_RCE[n_points,*])

	result_roi  = Estimate_RCE_HoffmannModel(ARR_RCE=arr_RCE_mean[nframe_injection:nframe_end],$
		ARR_TIME=arr_time_adjust, ST_DATA=st_data_hoffmann, PCERROR=pcerror, $
		POSITIVE_VALUES=positive_values,QUIET=quiet, METHOD=opt_method, DATA_ERR=data_err, YFIT=yfit, NAMES_PARAM=names_param_t)

	IF MIN(result_roi[0:5]) LT 0 THEN BEGIN
		error_str = 'ROI estimation failed, please check ROI selection or initial values'
		error_id  = 10l
	ENDIF

	st_data_r = {kep_ini:result_roi[0],$
		kel_ini:result_roi[1],$
		ah_ini :result_roi[2],$
		tau_ini:result_roi[4],$
		Bh_ini :result_roi[5]}

	kep_est = result_roi[0] &	kel_est = result_roi[1]
	A_est   = result_roi[2] &	Akep_est= result_roi[3]
	tau_est = result_roi[4] &	B_est   = result_roi[5]
	;----------------------------------------------------------------------------------
	IF opt_method EQ 5 THEN BEGIN ; sum  AIF
		yfit+=(B_est*arr_aif)
	ENDIF
	;----------------------------------------------------------------------------------

	st_data_ini = {$
		kep_ini  : result_roi[0],$
		kel_ini  : result_roi[1],$
		Ah_ini   : result_roi[2],$
		tau_ini  : result_roi[4],$
		Akep_ini : result_roi[0]*result_roi[2],$
		Bh_ini   : result_roi[5],$
		tau_infusion:st_data.tau_infusion[0]}

	;--------------------------------------------------
	data_adjust[n_points,nframe_injection:nframe_end] = yfit
	im_param[n_points,*]  = result_roi
	im_pcerror[n_points,*]    = pcerror
	;--------------------------------------------------
	im_statistics[n_points,0:1] = (Calculate_Statistics(arr_RCE_mean[nframe_injection:nframe_end], yfit, STRARR_NAMES=strarr_names1))[0:1]
	im_statistics[n_points,2:3] = (Calculate_Statistics(arr_RCE_mean[nframe_injection:nframe_end], yfit, RANGE=stat_range, STRARR_NAMES=strarr_names2))[0:1]
	;--------------------------------------------------

	names_stat  = [strarr_names1[0:1], strarr_names2[0:1]+'(Range)']
	names_param = names_param_t


ENDIF ELSE BEGIN
	; Sin estimación previa de ROI
	st_data_ini = st_data_hoffmann
ENDELSE
;------------------------------------------------------------------------------------

str_n_points = STRTRIM(STRING(n_points),2)
FOR n=0l, n_points-1 DO BEGIN

	IF (n+1) MOD 100 EQ 0 THEN PRINT, STRTRIM(n+1,2) + ' of ' + str_n_points

	arr_RCE = (REFORM(data_signal_RCE[n,*]))[nframe_injection:nframe_end]
	;-------------------------------------------------------------------------------------------
	result  = Estimate_RCE_HoffmannModel(ARR_RCE=arr_RCE, ARR_TIME=arr_time_adjust, $
		DATA_ERR=data_err,ST_DATA=st_data_ini, $
		PCERROR=pcerror, QUIET=quiet,POSITIVE_VALUES=positive_values, $
		YFIT=yfit, BESTNORM=bestnorm, METHOD=opt_method)
	;-------------------------------------------------------------------------------------------
	data_adjust[n,nframe_injection:nframe_end] = yfit
	im_param[n,*]  = result
	im_pcerror[n,*]    = pcerror ; puede ser kep, kel, y Akep con method = 2 (y tau)

	;--------------------------------------------------
	im_statistics[n,0:1] = (Calculate_Statistics(arr_RCE, yfit))[0:1]
	im_statistics[n,2:3] = (Calculate_Statistics(arr_RCE, yfit, RANGE=stat_range))[0:1]
	;--------------------------------------------------

ENDFOR

RETURN, 1

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Estimate_RCE_Array_LarssonModel, data,$
	ARR_TIME =arr_time,   $
	DATA_ERR=data_err,    $
	ST_DATA=st_data,      $
	ST_CP_MODEL=st_cp_model, $
	IM_PARAM=im_param,    $
	IM_PCERRORS=im_pcerrors,     $
	IM_STATISTICS=im_statistics, $
	STAT_RANGE= stat_range,$
	NO_ROI_ESTIMATION=no_roi_estimation,$
	DATA_SIGNAL_RCE=data_signal_rce,       $
	DATA_ADJUST=data_adjust, $
	METHOD=method,$
	QUIET=quiet,  $
	NAMES_PARAM=names_param,$
	NAMES_STAT =names_stat,$
	ERROR_ID=error_id, ERROR_STR=error_str

;frame_injection - Frame of injection ..

;st_data = {a1,a2,m1,m2}

st_data_larsson = { $
	kep_ini:st_data.kepl_ini,$
	s_ini:st_data.sl_ini}


arr_aux = [st_cp_model.a1,st_cp_model.a2,st_cp_model.m1,st_cp_model.m2]

frame_time       = st_data.frame_period          ; minutos
nframe_injection = st_data.nframe_injection      ; donde empieza la inyección del bolo
nframes_adjust   = st_data.nframes_adjust

nframe_end     = (nframe_injection+nframes_adjust)

; data_signal_RCE son valores de contraste relativo en tanto por uno (ya deben de estar normalizados)
; El primer valor no tiene porqué ser uno.. (está en media con respecto a antes del contraste)

error_id = 0l & error_str = ''

IF KEYWORD_SET(method) THEN opt_method=FIX(method[0]) ELSE  opt_method = 1l

; Data_err es un valor estimado de DESVIACIÓN ESTÁNDAR de los datos..( para entrada ERR de mpfitfun)
; estima el modelo de Hoffmann para varios puntos
; data_signal_RCE en % por uno

n_parameters = 3l  ; {si, kep, tau_shift}
n_statistics = 4l ; RMSE,ME,  RMSE(range),ME (range)

CASE SIZE(REFORM(data), /N_DIMENSIONS) OF
	2 : BEGIN; Array [n_points, n_frames] (100,nframes)
 		type_data   = 'ARRAY'
 		n_points    = (SIZE(REFORM(data), /DIMENSIONS))[0]
 		nframes     = (SIZE(REFORM(data), /DIMENSIONS))[1]
 		im_param    = FLTARR(n_points+1, n_parameters)
 		im_pcerrors = FLTARR(n_points+1, n_parameters)
 		im_statistics   = FLTARR(n_points+1, n_statistics)
		Data_signal_RCE = FLTARR(n_points+1,nframes)
		data_adjust = FLTARR(n_points+1,nframes)
		data_adjust[*,0:nframe_injection]=1.0
	END
 	ELSE : BEGIN
 		error_id = 1 &	error_str = 'Not correct dimensions of data'
 		RETURN, -1
 	END
ENDCASE

arr_time = INDGEN(nframes)*1.0*frame_time
arr_time_adjust = arr_time[nframe_injection:nframe_end]-(frame_time*nframe_injection)


Data_signal_RCE[0:n_points-1,*]  = RCE_signalRel_from_Signal(data, FRAME_INJECTION=nframe_injection) ; only calculates RCE in % per one
Data_signal_RCE[n_points,*]      = RCE_signalRel_from_Signal(TOTAL(data,1)/n_points, FRAME_INJECTION=nframe_injection) ; mean value

;------------------------------------------------------------------------------------
IF KEYWORD_SET(no_roi_estimation) EQ 0 THEN BEGIN
	; Primer paso, estima los parámetros medios para todo el ROI

	arr_RCE_mean   = Data_signal_RCE[n_points,*]

	result_roi  = Estimate_RCE_LarssonModel(ARR_RCE=arr_RCE_mean[nframe_injection:nframe_end],$
		ARR_TIME=arr_time_adjust, PCERROR=pcerror, ST_DATA=st_data_larsson, ARR_AUX=arr_aux, ST_NAMES_PARAM=names_param_t,$
		QUIET=quiet, METHOD=opt_method, DATA_ERR=data_err, YFIT=yfit)

	IF MIN(result_roi[0:1]) LE 0 THEN BEGIN
		error_str = 'ROI estimation failed, please check ROI selection or initial values'
		error_id  = 10l
	ENDIF

	kep_est = result_roi[0] &	S_est = result_roi[1]

	st_data_ini = {$
		kep_ini  : result_roi[0],$
		S_ini    : result_roi[1],$
		a1  :st_data_larsson.a1,$
		a2:st_data_larsson.a2,$
		m1:st_data_larsson.m1,$
		m2:st_data_larsson.m2}
	;--------------------------------------------------
	im_param[n_points,0:1]  = result_roi
	data_adjust[n_points,nframe_injection:nframe_end] = yfit
	im_pcerrors[n_points,0:1]   =  pcerror ; puede ser kep, kel, y Akep con method = 2 (y tau)
	;--------------------------------------------------
	im_statistics[n_points,0:1] = (Calculate_Statistics(arr_RCE_mean[nframe_injection:nframe_end], yfit, STRARR_NAMES=strarr_names1))[0:1]
	im_statistics[n_points,2:3] = (Calculate_Statistics(arr_RCE_mean[nframe_injection:nframe_end], yfit, RANGE=stat_range, STRARR_NAMES=strarr_names2))[0:1]

	names_stat  = [strarr_names1[0:1], strarr_names2[0:1]+'(Range)']
	names_param = names_param_t
	;--------------------------------------------------

ENDIF ELSE BEGIN
	; Sin estimación previa de ROI
	st_data_ini = st_data_larsson
ENDELSE
;------------------------------------------------------------------------------------

str_n_points = STRTRIM(STRING(n_points),2)
FOR n=0l, n_points-1 DO BEGIN

	IF (n+1) MOD 100 EQ 0 THEN PRINT, STRTRIM(n+1,2) + ' of ' + str_n_points

	arr_RCE = (REFORM(data_signal_RCE[n,*]))[nframe_injection:nframe_end]
	;--------------------------------------------------------------------------------
	result  = Estimate_RCE_LarssonModel(ARR_RCE=arr_RCE, ARR_TIME=arr_time_adjust, $
		DATA_ERR=data_err,ST_DATA=st_data_ini, PCERROR=pcerror, QUIET=quiet, ARR_AUX=arr_aux, $
		YFIT=yfit, BESTNORM=bestnorm, METHOD=opt_method )
	;--------------------------------------------------------------------------------
	im_param[n,0:1]  = result
	data_adjust[n,nframe_injection:nframe_end] = yfit
	im_pcerrors[n,0:1]   =  pcerror ; puede ser kep, kel, y Akep con method = 2 (y tau)

	;--------------------------------------------------
	im_statistics[n,0:1] = (Calculate_Statistics(arr_RCE, yfit))[0:1]
	im_statistics[n,2:3] = (Calculate_Statistics(arr_RCE, yfit, RANGE=stat_range))[0:1]
	;--------------------------------------------------
ENDFOR

RETURN, 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Estimate_RCE_HoffmannModel, ARR_RCE=arr_RCE,$
	ARR_TIME=arr_time, $
	ST_DATA=st_data,   $
	DATA_ERR=data_err, $
	QUIET= quiet,$
	METHOD=method,$
	YFIT=yfit, $
	POSITIVE_VALUES=positive_values,$
	BESTNORM=bestnorm,$
	PCERROR=pcerror,$
	NAMES_PARAM=names_param

; Notación: (Tofts, 1997)

;Result: 6 parameters: [kep, kel, Ah, Akep, tau_shift, Bh]

; Methods:
; 1 - Hoffmann classic (kep, kel, Ah)
; 2 - With Akep  (kep, kel, Akep)
; 3 - With Tau shift
; 4 - With Tau infusion (Brix model)
; 5 - Término adicional (Bh) pero sin tau_shift ni tau infusion

; YFIT     - Curva ajustada
; BESTNORM - Suma de diferencias al cuadrado entre el modelo y los datos
; DATA_ERR - Desviación estándar de los datos
;            (Asumiremos que es igual en todas las muestras, por eso se proporciona solo un valor)

mpside_opt=2 ; por defecto 0.

opt_positive = KEYWORD_SET(positive_values)

IF KEYWORD_SET(method) THEN opt_method=FIX(method[0]) ELSE  opt_method = 1l

n_frames = N_ELEMENTS(arr_time)
IF N_ELEMENTS(arr_RCE) NE n_frames THEN RETURN, -1

IF N_ELEMENTS(data_err) NE 0 THEN $
	rerr = REPLICATE(data_err[0], n_frames)

kep_ini = st_data.kep_ini[0]
kel_ini = st_data.kel_ini[0]
Ah_ini   = st_data.Ah_ini[0]

result_total = FLTARR(6)   ; 6  valores: [kep, kel, Ah, Akep, tau, Bh]
perror_total = FLTARR(6)   ; 6  valores de l-sigma error: [kep, kel, Ah, Akep, tau, Bh] (0, cuando no se calcula)

names_param = ['KEP', 'KEL', 'AH', 'AKEP', 'TAU_SHIFT', 'BH']

CASE opt_method OF
1 : BEGIN ; estándar, estimación de tres parámetros (kep, kel, A)
	X = FLOAT(arr_time)
	R = FLOAT(arr_RCE)
	start = FLOAT([kep_ini, kel_ini, Ah_ini])

	IF opt_positive THEN BEGIN
		parinfo =REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0], mpside:mpside_opt}, 3)
		parinfo[*].value = start ;[kep, kel, A]
		; los tres parametros limitados cantidades positivas
	ENDIF
	result = MPFITFUN('MPFIT_RCE_HOFFMANN_V1', X, R, rerr, start, QUIET=quiet, $
		PARINFO=parinfo, YFIT=yfit, COVAR=covar, PERROR=perror, BESTNORM=bestnorm)
	;kep = result[0];kel = result[1]	;A  = result[2]
	result_total[0:2] = result[0:2]
	result_total[3]   = result[2]*result[0]

	IF N_ELEMENTS(perror) EQ 3 THEN BEGIN
		perror_total[0:2] = perror[0:2]
		; Akep = a*kep
		perror_total[3] = StandardError_multiply(MEANS=DOUBLE([result[0],result[2]]), SD=DOUBLE([perror[0], perror[2]]))
	ENDIF

END
2 : BEGIN ; estimación directa del parámetro Akep (kep, kel, Akep)
	X = FLOAT(arr_time)
	R = FLOAT(arr_RCE)
	start = FLOAT([kep_ini, kel_ini, Ah_ini*kep_ini])

	IF opt_positive THEN BEGIN
		parinfo =REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0], mpside:mpside_opt}, 3)
		parinfo[*].value = start ;[kep, kel, Akep]
		; los tres parametros limitados cantidades positivas
	ENDIF
	result = MPFITFUN('MPFIT_RCE_HOFFMANN_V2', X, R, rerr, start, QUIET=quiet,$
		PARINFO=parinfo, YFIT=yfit, COVAR=covar, PERROR=perror, BESTNORM=bestnorm)
	;kep = result[0]	;kel = result[1]	;Akep   = result[2]
	result_total[0:1] = result[0:1]
	result_total[2]   = result[2]/result[0]; extracts A from A*kep
	result_total[3]   = result[2]
	IF N_ELEMENTS(perror) EQ 3 THEN BEGIN
		perror_total[0:1] = perror[0:1]
		perror_total[3]   = perror[2]
		perror_total[2]   = StandardError_divide(MEANS=DOUBLE([result[2],result[0]]), SD=DOUBLE([perror[2], perror[0]]))
	ENDIF

END
3 : BEGIN ; Estimación de cuatro parámetros (kep, kel, A, tau_shift)

	tau_ini = st_data.tau_ini

	X = FLOAT(arr_time)
	R = FLOAT(arr_RCE)
	start = FLOAT([kep_ini, kel_ini, Ah_ini, tau_ini])

	IF opt_positive THEN BEGIN
		parinfo =REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0], mpside:mpside_opt}, 4)
		parinfo[*].value = start ;[kep, kel, A, tau]
		; los tres parametros limitados cantidades positivas
		parinfo[3].limited = [0,0] ; don't limit the time

	ENDIF
	result = MPFITFUN('MPFIT_RCE_HOFFMANN_V3', X, R, rerr, start, QUIET=quiet,$
		PARINFO=parinfo, YFIT=yfit, COVAR=covar, PERROR=perror, BESTNORM=bestnorm)
	;kep = result[0]	;kel = result[1]	;A   = result[2]	;tau_inf = result[3]
	result_total[0:2] = result[0:2]
	result_total[3]   = result[2]*result[0]; extracts Akep from A and kep
	result_total[4]   = result[3]
	IF N_ELEMENTS(perror) EQ 4 THEN BEGIN
		perror_total[0:2] = perror[0:2]
		perror_total[3]   = StandardError_multiply(MEANS=DOUBLE([result[2],result[0]]), SD=DOUBLE([perror[2], perror[0]]))
		perror_total[4]   = perror[3]
	ENDIF

END
4 : BEGIN
	; tau_infusion; tiepo de infusión (pequeño ...)

	tau_infusion = st_data.tau_infusion[0]

	X = FLOAT([tau_infusion[0],arr_time])
	R = FLOAT(arr_RCE)
	start = FLOAT([kep_ini, kel_ini, Ah_ini])
	IF opt_positive THEN BEGIN
		parinfo =REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0], mpside:mpside_opt}, 3)
		parinfo[*].value = start ;[kep, kel, A]
	; los tres parametros limitados cantidades positivas
	ENDIF
	result = MPFITFUN('MPFIT_RCE_HOFFMANN_V4', X, R, rerr, start, QUIET=quiet,$
		PARINFO=parinfo, YFIT=yfit, COVAR=covar, PERROR=perror, BESTNORM=bestnorm)
	;result =; kep, kel, A
	result_total[0:2] = result[0:2]
	result_total[3]   = result[2]*result[0]
	IF N_ELEMENTS(perror) EQ 3 THEN BEGIN
		perror_total[0:2] = perror[0:2]
		perror_total[3]   = StandardError_multiply(MEANS=DOUBLE([result[2],result[0]]), SD=DOUBLE([perror[2], perror[0]]))
	ENDIF
END
5 : BEGIN
	; prueba adiciona de incluir una parte de suma de la AIF (contribución de voxeles vascularizados)
	Bh_ini = st_data.Bh_ini ;

	X = FLOAT([arr_time, arr_AIF])
	R = FLOAT(arr_RCE)
	start = FLOAT([kep_ini, kel_ini, Ah_ini, Bh_ini])
	IF opt_positive THEN BEGIN
		parinfo =REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0], mpside:mpside_opt}, 4)
		parinfo[*].value = start ;[kep, kel, A, B]
	; los tres parametros limitados cantidades positivas
	ENDIF
	result = MPFITFUN('MPFIT_RCE_HOFFMANN_V5', X, R, rerr, start, QUIET=quiet,$
		PARINFO=parinfo, YFIT=yfit, COVAR=covar, PERROR=perror, BESTNORM=bestnorm)
	;result =; kep, kel, A
	result_total[0:2] = result[0:2]
	result_total[3]   = result[2]*result[0]
	result_total[5]   = result[3]
	IF N_ELEMENTS(perror) EQ 4 THEN BEGIN
		perror_total[0:2] = perror[0:2]
		perror_total[3]   = StandardError_multiply(MEANS=DOUBLE([result[2],result[0]]), SD=DOUBLE([perror[2], perror[0]]))
		perror_total[5]   = perror[3]
	ENDIF
END
ELSE : RETURN, -1
ENDCASE

degrees_of_freedom = N_ELEMENTS(R)-N_ELEMENTS(start)
pcerror = perror_total*SQRT(bestnorm/degrees_of_freedom)

;-------------------------------------------
; Infinite or not defined values
pos = WHERE(FINITE(result) NE 1,ct)
IF ct GT 0 THEN result[pos]=-1
;-------------------------------------------

RETURN, result_total

END


;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Estimate_RCE_LarssonModel, ARR_RCE=arr_RCE, $
	ARR_TIME=arr_time, $
	ST_DATA =st_data,  $
	ARR_AUX =arr_aux,  $
	DATA_ERR=data_err, $
	QUIET=quiet,       $
	METHOD=method,     $
	YFIT=yfit,         $
	BESTNORM=bestnorm, $
	PCERROR=pcerror,   $
	NAMES_PARAM=names_param

; Ajuste de (Larsson et al, 1990), fórmula descrita en (Tofts et al, 1997) (ecuación 16)
;

;arr_ini = [kep, Si] ;
;Result: 2 parameters: [kep, Si]

; YFIT     - Curva ajustada
; BESTNORM - Suma de diferencias al cuadrado entre el modelo y los datos
; DATA_ERR - Desviación estándar de los datos
;            (Asumiremos que es igual en todas las muestras, por eso se proporciona solo un valor)


IF KEYWORD_SET(method) THEN opt_method=FIX(method[0]) ELSE  opt_method = 1l

n_frames = N_ELEMENTS(arr_time)
IF N_ELEMENTS(arr_RCE) NE n_frames THEN RETURN, -1

IF N_ELEMENTS(data_err) NE 0 THEN $
	rerr = REPLICATE(data_err[0], n_frames)

result_total = FLTARR(2)   ; 2  valores: [kep, Si]
perror_total = FLTARR(2)   ; 2  valores de l-sigma error:  (0, cuando no se calcula)

perror_rce = FLTARR(3)

kep_ini = st_data.kep_ini[0]
S_ini   = st_data.S_ini[0]

names_param = ['KEP', 'SL']

a1 = arr_aux[0]
a2 = arr_aux[1]
m1 = arr_aux[2]
m2 = arr_aux[3]

mpside_opt = 2; por defecto, 0

	X = FLOAT([a1,a2,m1,m2,arr_time])
	R = FLOAT(arr_RCE)
	start = FLOAT([kep_ini, S_ini])

	parinfo = REPLICATE({value:0.D, fixed:0, limited:[1,0], limits:[0.D,0], mpside:mpside_opt}, 2)
	parinfo[*].value = FLOAT([kep_ini, S_ini]) ;
	parinfo[0:1].limited = [0,1] ; limitado a maximo uno (valor negativo)
	; los tres parametros limitados cantidades positivas

	;quiet = 0
	result = MPFITFUN('MPFIT_RCE_LARSSON_V1', X, R, rerr, start, QUIET=quiet, $;PARINFO=parinfo, $
		YFIT=yfit, COVAR=covar, PERROR=perror, BESTNORM=bestnorm)
	;kep = result[0];kel = result[1]	;A  = result[2]
	kep_est = result[0]
	Si_est  = result[1]

	result_total=[kep_est, Si_est]

	IF N_ELEMENTS(perror) EQ N_ELEMENTS(result) THEN BEGIN
		perror_total[0:1] = perror[0:1]
	ENDIF


degrees_of_freedom = N_ELEMENTS(R)-N_ELEMENTS(start)
pcerror = perror_total*SQRT(bestnorm/degrees_of_freedom)

;-------------------------------------------
; Infinite or not defined values
pos = WHERE(FINITE(result) NE 1,ct)
IF ct GT 0 THEN result[pos]=-1
;-------------------------------------------

RETURN, result_total

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Filter_Gaussian_1D, array, SIGMA=sigma, FLOATING=floating

;; Function filter_Gaussian_1d.pro
;;
;; Resultado, LONG o FLOAT
;;

IF N_PARAMS() NE 1 THEN BEGIN & PRINT, 'Datos no válidos'    & RETURN, -1 & END

IF NOT KEYWORD_SET(floating) THEN but_float = 0 ELSE but_float = 1
IF N_ELEMENTS(sigma) EQ 0 THEN sigma_f = 1.0  ELSE sigma_f = FLOAT(sigma[0])

;IF (SIZE(array))[0] NE 1    THEN BEGIN & PRINT, 'Datos no válidos'    & RETURN, -1 & END

array_o = FLOAT(array)

sigma_f = sigma_f > 0.5d

flt_arr = INDGEN(LONG((sigma+1)^2) + ((LONG((sigma+1)^2) MOD 2) EQ 0))
tam     = N_ELEMENTS(flt_arr)
flt_arr = (FLOAT(flt_arr - (tam/2)))^2
struct = EXP(-flt_arr/(2*(sigma_f^2)))

struct = struct[WHERE(struct GT 0.05)]
tam = N_ELEMENTS(struct)

struct1 = REFORM(struct, tam)
scl     = TOTAL(struct)
array_convol = CONVOL(array_o, struct1, scl, /CENTER, /EDGE_WRAP)

IF but_float EQ 1 THEN BEGIN
	RETURN, FLOAT(array_convol)
ENDIF ELSE BEGIN
	RETURN, LONG(ROUND(array_convol))
ENDELSE

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION finite_values, array, MAX_VALUE=max_value, MIN_VALUE=min_value

IF N_ELEMENTS(max_value) EQ 0 THEN max_value = (MACHAR()).xmax/2
IF N_ELEMENTS(min_value) EQ 0 THEN min_value = (MACHAR()).xmin/2

n_elem = N_ELEMENTS(array)

arr_sal = array

pos = WHERE(FINITE(arr_sal, /INF, SIGN=1) EQ 1, ct)
IF ct GE 1 THEN BEGIN
	arr_sal[pos] = max_value
ENDIF
pos = WHERE(FINITE(arr_sal, /INF, SIGN=-1) EQ 1, ct)
IF ct GE 1 THEN BEGIN
	arr_sal[pos] = min_value
ENDIF

pos = WHERE(FINITE(arr_sal, /NAN) EQ 1, ct)
IF ct GE 1 THEN BEGIN
	arr_sal[pos] = 0
ENDIF

RETURN, arr_sal

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<++>
;+
; NAME:
;       FSC_NORMALIZE
;
; PURPOSE:
;
;       This is a utility routine to calculate the scaling vector
;       required to position a graphics primitive of specified range
;       at a specific position in an arbitray coordinate system. The
;       scaling vector is given as a two-element array like this:
;
;          scalingVector = [translationFactor, scalingFactor]
;
;       The scaling vector should be used with the [XYZ]COORD_CONV
;       keywords of a graphics object or model. For example, if you
;       wanted to scale an X axis into the coordinate range of -0.5 to 0.5,
;       you might type something like this:
;
;          xAxis->GetProperty, Range=xRange
;          xScale = FSC_Normalize(xRange, Position=[-0.5, 0.5])
;          xAxis, XCoord_Conv=xScale
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: david@idlcoyote.com
;       Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; CATEGORY:

;       Object Graphics
;
; CALLING SEQUENCE:
;
;       xscaling = FSC_NORMALIZE(xrange, POSITION=position)
;
; INPUTS:
;
;       XRANGE: A two-element vector specifying the data range.
;
; KEYWORD PARAMETERS:
;
;       POSITION: A two-element vector specifying the location
;       in the coordinate system you are scaling into. The vector [0,1]
;       is used by default if POSITION is not specified.
;
; COMMON BLOCKS:
;
;       None.
;
; EXAMPLE:
;
;       See above.
;
; MODIFICATION HISTORY:
;       Written by:  David W. Fanning, OCT 1997.
;       Fixed a problem with illegal divide by zero. 21 April 2005. DWF.
;       Fixed a problem when range[0] is greater than range[1]. 11 July 2006. DWF.
;       Renamed to FSC_Normalize to avoid conflicts with 10,000 other programs named NORMALIZE. 17 October 2008. DWF.
;-
;******************************************************************************************;
;  Copyright (c) 2008, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;

FUNCTION FSC_Normalize, range, Position=position

On_Error, 1
IF N_Params() EQ 0 THEN Message, 'Please pass range vector as argument.'

IF (N_Elements(position) EQ 0) THEN position = [0.0D, 1.0D] ELSE $
    position=Double(position)
range = Double(range)

IF range[1] GE range[0] THEN BEGIN
   scale = [((position[0]*range[1])-(position[1]*range[0])) / $
       ((range[1]-range[0]) > 1e-12), (position[1]-position[0])/((range[1]-range[0]) > 1e-12)]
ENDIF ELSE BEGIN
   scale = [((position[1]*range[0])-(position[0]*range[1])) / $
       ((range[0]-range[1]) > 1e-12), (position[1]-position[0])/((range[0]-range[1]) > 1e-12)]
   scale[1] = -scale[1]
ENDELSE
RETURN, scale
END
;-------------------------------------------------------------------------
;<-->
;<+>

PRO  fwidget_ControlZoomDraw, widget_id_draw, REL_XY=rel_xy, ZOOM=zoom, XLIMITS=xlimits, YLIMITS=ylimits, SCROLL_LIMIT=scroll_limit

	IF N_ELEMENTS(rel_xy) EQ 0 THEN rel_xy = 1d ; aspect relation between image size in X and Y

	zoomd   = DOUBLE(zoom)
	st_draw = WIDGET_INFO(widget_id_draw, /GEOMETRY)

	WIDGET_CONTROL,  widget_id_draw, GET_VALUE = owindow
	owindow->getproperty, CURRENT_ZOOM=current_zoom, LOCATION=location, VISIBLE_LOCATION=visible_location
	PRINT, visible_location

	draw_xsize =st_draw.draw_xsize
	draw_ysize =st_draw.draw_ysize
	xsize = st_draw.xsize
	ysize = st_draw.ysize

	IF ((st_draw.draw_xsize*zoomd GE scroll_limit) OR (st_draw.draw_ysize*zoomd GE scroll_limit)) THEN RETURN ; absolute limit cannot be exceeded

	IF (zoomd GT 1) THEN BEGIN
		lt_x_max = st_draw.xsize*zoomd LT xlimits[1]
		lt_y_max = st_draw.ysize*zoomd LT ylimits[1]
		IF (lt_x_max EQ 1) THEN xsize = st_draw.xsize*zoomd   ; values can be separated, because the deformation is due to draw_xsize and draw_ysize
		IF (lt_y_max EQ 1) THEN ysize = st_draw.ysize*zoomd

		lt_x_maxS = st_draw.draw_xsize*zoomd LT scroll_limit[0]
		lt_y_maxS = st_draw.draw_ysize*zoomd LT scroll_limit[0]

		IF (lt_x_maxS EQ 1) AND (lt_y_maxS EQ 1) THEN BEGIN
			draw_xsize = st_draw.draw_xsize*zoomd ; scroll size must be always be joined or image will be deformed
			draw_ysize = st_draw.draw_ysize*zoomd
		ENDIF
		;WIDGET_CONTROL, widget_id_draw, XSIZE=xsize, YSIZE=ysize, DRAW_YSIZE=draw_ysize, DRAW_XSIZE=draw_xsize

	ENDIF ELSE BEGIN ; reducing image

		gt_x_minS = st_draw.draw_xsize*zoomd  GT xlimits[0]
		gt_y_minS = st_draw.Draw_ysize*zoomd  GT ylimits[0]

		IF gt_x_minS AND gt_y_minS THEN BEGIN
			draw_xsize = st_draw.draw_xsize*zoomd
			draw_ysize = st_draw.draw_ysize*zoomd  ; values cannot be separated, and scroll option is reduced before
		ENDIF
		gt_x_min = st_draw.xsize*zoomd GT xlimits[0]
		gt_y_min = st_draw.ysize*zoomd GT ylimits[0]

		IF gt_x_min THEN BEGIN
			IF st_draw.xsize GT draw_xsize THEN xsize = draw_xsize
		ENDIF
		IF gt_y_min THEN BEGIN
			IF st_draw.ysize GT draw_ysize THEN ysize = draw_ysize
		ENDIF
		;WIDGET_CONTROL, widget_id_draw, XSIZE=xsize, YSIZE=ysize, DRAW_YSIZE=draw_ysize, DRAW_XSIZE=draw_xsize
	ENDELSE

	WIDGET_CONTROL, widget_id_draw, XSIZE=xsize, YSIZE=ysize, DRAW_YSIZE=draw_ysize, DRAW_XSIZE=draw_xsize

	;owindow->setproperty, VISIBLE_LOCATION=[xsize/4, ysize/4]

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Get_BITUPM_strarray
; Retorna un array (en vertical, con el nombre del copyright)

strarr_info = STRARR(1,3)

strarr_info[0] = 'BIT (Biomedical image Technologies)'
strarr_info[1] = 'Departamento de ingeniería Electrónica'
strarr_info[2] =  STRING(169B) + ' Universidad Politécnica de Madrid'

RETURN, strarr_info

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Get_color, color, HEX=hex, BYTE=byte

; Traduce colores de hexadecimal a 3 valores byte y viceversa

type = SIZE(color, /TYPE)
IF type EQ 3 THEN BEGIN
    RETURN, STRTRIM(STRING(color, FORMAT='(z)'),2)
ENDIF
IF N_ELEMENTS(color) EQ 3 THEN BEGIN ; En formato de tres números
   col = LONG(color[0]) + LONG(color[1])*256l + LONG(color[2])*256l*256l
ENDIF

RETURN, col
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION  Get_complete_path, path

IF N_ELEMENTS(path) EQ 0  THEN RETURN, ''
IF SIZE(path, /TYPE) ne 7 THEN RETURN, -1

; Añande un "\" a un STRING (que representa un path) si le falta al final.

n_files = N_ELEMENTS(path)
IF n_files EQ 1 THEN path_mark = STRARR(1) ELSE path_mark = STRARR(n_files)
FOR i=0l, n_files-1 DO BEGIN
	path_mark[i] = FILE_DIRNAME(path[i] + PATH_SEP() + 'A', /MARK_DIRECTORY)
ENDFOR
IF n_files EQ 1 THEN path_mark = path_mark[0]
RETURN, path_mark

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Get_disc, diameter, QUADRANT=quadrant, OCTANT=octant, DIAGONAL=diagonal

; Crea un disco de unos (elemneto estructurante típico) de diametro determinado
; Por lo tanto, el disco tendrá diámetro PAR o IMPAR

image = BYTARR(diameter, diameter)
xsize = DOUBLE(diameter)

rad   = DOUBLE(diameter)/2d

FOR x=0, xsize-1 DO BEGIN
    FOR y=0, xsize-1 DO BEGIN
        IF (((x-rad+0.5)^2 + (y-rad+0.5)^2) LT ((rad-0.1)^2)) THEN $
        image[x,y] = 1
    ENDFOR
ENDFOR

IF KEYWORD_SET(quadrant) THEN BEGIN
	image[*,diameter/2:diameter-1] = 0
	image[diameter/2:diameter-1,*] = 0
ENDIF

IF KEYWORD_SET(octant) THEN BEGIN
	image[*,diameter/2:diameter-1] = 0
	image[diameter/2:diameter-1,*] = 0
	p = WHERE(image GT 0, ct)
	p_x = p MOD diameter
	p_y = p / diameter
	pos = WHERE(p_y GT p_x, ct2)
	image[p[pos]]=0
ENDIF

IF KEYWORD_SET(diagonal) THEN BEGIN
	p = WHERE(image GT 0, ct)
	p_x = p MOD diameter
	p_y = p / diameter
	pos = WHERE(p_y NE p_x, ct)
	image[p[pos]]=0
ENDIF

RETURN, image

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Get_elapsedtime, time_initial
; Retorna el tiempo transcurrido desde la entrada ,time_initial (valor en segundos)

time_actual = SYSTIME(1)
time_elapsed = time_actual - time_initial
RETURN, time_elapsed
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Get_FloatingFromString, str, NEW_STRING=new_string, ERROR_ID=error_id, $
	POSITIVE=positive, MAX_DECIMALS=max_decimals

error_id = 0l

opt_positive  = KEYWORD_SET(positive)

len_str = STRLEN(str)
IF len_str EQ 0 THEN BEGIN
	new_string = ''
	IF (opt_positive EQ 0) THEN BEGIN
		error_id = 1 & RETURN, -1
	ENDIF
ENDIF

first_char = BYTE(STRMID(STRTRIM(str,2),0,1))
IF first_char EQ 45 THEN BEGIN ;45 = '-'
	IF (opt_positive EQ 0) THEN BEGIN
		new_string = ''
	ENDIF ELSE BEGIN
		opt_negative = 1l
		first_char = BYTE(STRMID(STRTRIM(str,2),1,1))
	ENDELSE
ENDIF
IF (first_char LT 48) OR (first_char GT 57) THEN BEGIN
	IF (opt_positive EQ 0) AND (opt_posorzero EQ 0) THEN BEGIN
		error_id = 2 & new_string = '' & RETURN, -1
	ENDIF ELSE BEGIN
		new_string = ''
	ENDELSE
ENDIF

str_1 = STRMID(STRING(FIX(max_decimals[0])),2)
str_format = '(F20.'+ str_1 + ')'

number = FLOAT(str[0])

IF opt_positive   THEN number = number > 0

new_string = STRTRIM(STRING(number[0], FORMAT=str_format),2)
RETURN, number

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Get_IntegerFromString, str, NEW_STRING=new_string, ERROR_ID=error_id, $
	POSITIVE=positive, POSORZERO=posorzero, MAXVAL=maxval

error_id = 0l

opt_positive  = KEYWORD_SET(positive)
opt_posorzero = KEYWORD_SET(posorzero)

len_str = STRLEN(str)
IF len_str EQ 0 THEN BEGIN
	new_string = ''
	IF (opt_positive EQ 0) AND (opt_posorzero EQ 0) THEN BEGIN
		error_id = 1 & RETURN, -1
	ENDIF
ENDIF

first_char = BYTE(STRMID(STRTRIM(str,2),0,1))
IF first_char EQ 45 THEN BEGIN ;45 = '-'
	IF (opt_positive EQ 0) AND (opt_posorzero EQ 0) THEN BEGIN
		new_string = ''
	ENDIF ELSE BEGIN
		opt_negative = 1l
		first_char = BYTE(STRMID(STRTRIM(str,2),1,1))
	ENDELSE
ENDIF
IF (first_char LT 48) OR (first_char GT 57) THEN BEGIN
	IF (opt_positive EQ 0) AND (opt_posorzero EQ 0) THEN BEGIN
		error_id = 2 & new_string = '' & RETURN, -1
	ENDIF ELSE BEGIN
		new_string = ''
	ENDELSE
ENDIF
integer = LONG(str[0])

IF N_ELEMENTS(maxval) NE 0 THEN BEGIN
	IF integer GT  maxval[0] THEN integer = LONG(maxval[0])
	IF integer LT -maxval[0] THEN integer = LONG(-maxval[0])
ENDIF

IF opt_positive   THEN integer  = integer > 1
IF opt_posorzero  THEN integer = integer > 0


new_string = STRTRIM(STRING(integer[0]),2)

RETURN, integer

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Get_License_Array

Message_Text = [$
	'DCE@urLAB', $
'Copyright (C) 2012,2013, Universidad Politécnica de Madrid', $
'This program comes with ABSOLUTELY NO WARRANTY;', $
'This is free software, and you are welcome to redistribute it', $
'under certain conditions; read  "disclaimer.txt" file for details.']

Message_Text = TRANSPOSE(message_text)

RETURN, Message_text

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Get_name, file
; Función que extrae el nombre de fichero de un path + fichero

IF N_ELEMENTS(file) EQ 0  THEN RETURN, ''
IF SIZE(file, /TYPE) ne 7 THEN RETURN, -1

n_files = N_ELEMENTS(file)
str_s = PATH_SEP()
IF N_ELEMENTS(file) EQ 1 THEN files_sal = STRARR(1) ELSE files_sal = STRARR(SIZE(file, /DIMENSIONS))
FOR i=0l, n_files-1 DO BEGIN
	pos_sep = STRPOS(file[i], str_s, /REVERSE_SEARCH)
	IF pos_sep[0] EQ -1 THEN files_sal[i] = file[i] ELSE files_sal[i] = STRMID(file[i], pos_sep+1)
ENDFOR
IF N_ELEMENTS(files_sal) EQ 1 THEN files_sal = files_sal[0]
RETURN, files_sal
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Get_name_extension, file

; Función que extrae la extensión de un fichero (si no tiene, da '')

IF N_ELEMENTS(file) EQ 0  THEN RETURN, ''
IF SIZE(file, /TYPE) ne 7 THEN RETURN, -1


name = Get_name(file)
n_files = N_ELEMENTS(file)
str_p = '.'

arr_name_ext = STRARR(n_files)
FOR i=0l, n_files-1 DO BEGIN
	pos_sep = STRPOS(name[i], str_p, /REVERSE_SEARCH)
	IF pos_sep EQ -1 THEN BEGIN
		arr_name_ext[i] = ''
	ENDIF ELSE BEGIN
		name_extension = STRMID(name[i], pos_sep+1)
		IF N_ELEMENTS(name_extension) EQ 1 THEN arr_name_ext[i]= name_extension[0] ELSE BEGIN
			arr_name_ext[i] = name_extension[N_ELEMENTS(name_extension)-1]
		ENDELSE
	ENDELSE
ENDFOR
IF n_files EQ 1 THEN arr_name_ext = arr_name_ext[0]

RETURN, arr_name_ext
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Get_name_field, file

; Función que extrae el nombre de un fichero (sin extensión)

IF N_ELEMENTS(file) EQ 0  THEN RETURN, ''
IF SIZE(file, /TYPE) ne 7 THEN RETURN, -1


name = Get_name(file)
str_p = '.'
n_files = N_ELEMENTS(file)

arr_name_field = STRARR(n_files)

FOR i=0l, n_files-1 DO BEGIN
	pos_sep = STRPOS(name[i], str_p, /REVERSE_SEARCH)
	IF pos_sep EQ -1 THEN BEGIN
		arr_name_field[i] = name[i]
	ENDIF ELSE BEGIN
		name_field = STRMID(name[i], 0, pos_sep)
		IF N_ELEMENTS(name_field) EQ 1 THEN arr_name_field[i] = name_field[0]
	ENDELSE
ENDFOR
IF n_files EQ 1 THEN arr_name_field = arr_name_field[0]

RETURN, arr_name_field
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Get_OS_info

os      = !VERSION.OS
arch    = !VERSION.ARCH
release = !VERSION.RELEASE

RETURN, ['IDL ' + release + ' ' +os +  ' ' + arch]

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION  Get_path, file

IF N_ELEMENTS(file) EQ 0  THEN RETURN, ''
IF SIZE(file, /TYPE) ne 7 THEN RETURN, -1

path_temp = FILE_DIRNAME(file)
Path = Get_complete_path(path_temp)
RETURN, path

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION get_point_YdownConvention, point, DIMENSIONS=dimensions, SIZE_ROI=size_roi

;transforma un punto de convención Y-UP (esquina superior derecha) a convención IDL (esquina inferior derecha)

IF N_ELEMENTS(size_roi) EQ 0 THEN dim_roi = FIX([1,1]) ELSE dim_roi = size_roi[0:1]
pini = INTARR(2)
pini[0] = point[0]
pini[1] = dimensions[1]-point[1]-dim_roi[1]
RETURN ,pini

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION get_point_YupConvention, point, DIMENSIONS=dimensions, SIZE_ROI=size_roi

;transforma un punto de convención IDL (esquina inferior derecha) a convención Y-UP (esquina superior derecha)

IF N_ELEMENTS(size_roi) EQ 0 THEN dim_roi = FIX([1,1]) ELSE dim_roi = size_roi[0:1]
pini = INTARR(2)
pini[0] = point[0]
pini[1] = dimensions[1]-point[1]-dim_roi[1]
RETURN ,pini

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Get_ROIfromTxtdata, str_array, DIMENSIONS=dimensions, RESOLUTION=resolution,ERROR_ID=error_id

error_id=0l


IF N_ELEMENTS(dimensions) NE 2 THEN BEGIN
	str_dimensions = STRTRIM(((STRSPLIT(str_array[0], ':', /EXTRACT))[1]),2)
	str_dimensions = STRSPLIT(str_dimensions, ' ',/EXTRACT)
	xsize = FIX(str_dimensions[0])
	ysize = FIX(str_dimensions[1])
	dimensions = [xsize, ysize]
ENDIF ELSE BEGIN
	xsize = dimensions[0]
	ysize = dimensions[1]
ENDELSE

IF N_ELEMENTS(resolution) NE 1 THEN BEGIN
	str_resolution = STRTRIM(((STRSPLIT(str_array[1], ':', /EXTRACT))[1]),2)
	resolution     = FIX(str_resolution[0])
ENDIF


pos_comment = STRPOS(STRTRIM(str_array,2),'#')
pos_number  = WHERE(pos_comment EQ -1)

str_array = str_array[pos_number]

n_points = N_ELEMENTS(str_array)

opt_delete = 0l
IF opt_delete THEN BEGIN
	pos_out = WHERE(RANDOMU(seed,n_points) LT 0.8)
	str_array = str_array[pos_out]
	n_points = N_ELEMENTS(str_array)
ENDIF

arr_pos_xy = INTARR(2, n_points)
str_tab    = STRING(9B)

FOR i=0l, n_points-1 DO BEGIN
	arr_split = FIX(STRSPLIT(str_array[i], ' '+str_tab, /EXTRACT))
	IF N_ELEMENTS(arr_split) NE 2 THEN BEGIN
		error_id=2l
		RETURN, -1
	ENDIF
	arr_pos_xy[*,i] = arr_split
ENDFOR

;------------------------------------------------------------------------
; Pasa a convenio de coordenadas "natural" (desde 0, y con indice Y creciente)
arr_pos_xy-=1 ; resta
arr_pos_xy[1,*] = ysize-1-arr_pos_xy[1,*]


im_mask = BYTARR(xsize,ysize)

FOR i=0l, n_points-1 DO BEGIN
	x_ini = arr_pos_xy[0,i]
	y_end = arr_pos_xy[1,i]
	x_end = x_ini+resolution-1
	y_ini = y_end-resolution+1
	im_mask[x_ini:x_end,y_ini:y_end]=1
ENDFOR

RETURN, im_mask

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION get_windowNumber, OPTNEXT=optnext, ACTUAL=actual, XPOS=xpos, YPOS=ypos

ysize = 500

IF KEYWORD_SET(optnext) THEN BEGIN
	win_number=(!D.WINDOW+1) MOD 32

	DEVICE, GET_WINDOW_POSITION=wp, GET_SCREEN_SIZE = screenSize
	xpos = (screenSize[0]-20)  <(wp[0]-40) > 20
	ypos = (screenSize[1]-120) <(wp[1]+20-ysize) > 20

ENDIF ELSE BEGIN
	win_number=!D.WINDOW
	CASE win_number OF
	-1 : win_number = 1
	0  : BEGIN ; only main window...
		DEVICE,GET_WINDOW_POSITION=wp, GET_SCREEN_SIZE = screenSize
		xpos = (screenSize[0]-20)  <(wp[0]) > 20
		ypos = (screenSize[1]-120) <(wp[1]-ysize) > 20
	END
	ELSE: BEGIN
		DEVICE,GET_WINDOW_POSITION=wp, GET_SCREEN_SIZE = screenSize
		xpos = (screenSize[0]-20)  <(wp[0]-24)    > 20
		ypos = (screenSize[1]-120) <(wp[1]-ysize) > 20
    END
    ENDCASE
ENDELSE

RETURN, win_number
END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION GetBrukerHeader_Array, str_array, str_key, ERROR_ID=error_id

error_id = -1

n_lines = N_ELEMENTS(str_array)
FOR i=0, n_lines-1 DO BEGIN
	pos = STRPOS(str_array[i], str_key)
	IF pos EQ 0 THEN BREAK
ENDFOR
IF pos NE 0 THEN RETURN, -1

line = str_array[i]
line = STRSPLIT(line, '=', /EXTRACT)
IF N_ELEMENTS(line) NE 2 THEN RETURN, -2
line = STRTRIM(line[1],2)
line = STRSPLIT(line, '(', /EXTRACT, /PRESERVE_NULL)
IF N_ELEMENTS(line) NE 2 THEN RETURN, -3
line = STRTRIM(line[1],2)
line = STRSPLIT(line, ')', /EXTRACT, /PRESERVE_NULL)
IF N_ELEMENTS(line) NE 2 THEN RETURN, -4
line = STRTRIM(line[0],2)
number = FIX(line)
IF number LE 0 THEN RETURN, -2

n_numbers = 0l
arr_data = STRARR(number)

FOR j=1l, number DO BEGIN
	line = str_array[i+j]
	IF STRPOS(line, '##') NE -1 THEN BREAK
	arr_temp = STRSPLIT(line, ' ' , /EXTRACT)
	n_temp = N_ELEMENTS(arr_temp)
	arr_data[n_numbers:n_numbers+n_temp-1] = arr_temp
	n_numbers+=n_temp
	IF n_numbers GE number THEN BREAK
ENDFOR

IF n_numbers NE number THEN RETURN, -1

arr_data = FLOAT(arr_data)

error_id=0

RETURN, arr_data

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION GetBrukerHeader_key, str_array, str_key, ERROR_ID=error_id

error_id = -1

n_lines = N_ELEMENTS(str_array)
FOR i=0, n_lines-1 DO BEGIN
	pos = STRPOS(str_array[i], str_key)
	IF pos EQ 0 THEN BREAK
ENDFOR

IF pos NE 0 THEN RETURN, -1

line = str_array[i]
line = STRSPLIT(line, '=', /EXTRACT)
IF N_ELEMENTS(line) NE 2 THEN RETURN, -2
str_value = STRTRIM(line[1],2)


error_id = 0

RETURN, str_value

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION GetWidget_PositionInScreen, id_widget, CENTER=center, POSITION=position, PIXELS=pixels

; Coloca una widget en una posición EN LA PANTALLA
; Posibles posiciones:

;2       0;
;    c    ;  C= center
;3       1;

;    4    ;
;7   8   6;  C= center = 8
;    5    ;

IF N_ELEMENTS(id_widget) NE 1 THEN RETURN,  -1
IF SIZE(id_widget, /TYPE) NE 3 THEN RETURN, -2
IF WIDGET_INFO(id_widget, /VALID) NE 1 THEN RETURN, -3
IF KEYWORD_SET(center) THEN opt_center=1 ELSE opt_center=0
IF N_ELEMENTS(position) NE 0 THEN pos = FIX(position[0]) ELSE pos =8
IF opt_center THEN pos = 8

IF N_ELEMENTS(pixels) EQ 2 THEN pos = -1 ; los pixeles mandan...

DEVICE, Get_Screen_Size = screenSize
xCenter = screenSize[0]/2d
yCenter = screenSize[1]/2d

geom = WIDGET_INFO(id_widget, /GEOMETRY)
xsize = geom.Scr_XSize
ysize = geom.Scr_YSize
xHalfSize = geom.Scr_XSize/2d
yHalfSize = geom.Scr_YSize/2d

CASE pos OF
-1 : BEGIN & xoffset = pixels[0] & yoffset = pixels[1] & END
0: BEGIN & xoffset = screensize[0]-xsize  	& yoffset= 0  					& END
1: BEGIN & xoffset = screensize[0]-xsize  	& yoffset= screensize[1]-ysize  & END
2: BEGIN & xoffset = 0   					& yoffset= 0  					& END
3: BEGIN & xoffset = 0   					& yoffset= screensize[1]-ysize	& END
4: BEGIN & xoffset = xCenter-xHalfSize  	& yoffset= 0 					& END
5: BEGIN & xoffset = xCenter-xHalfSize   	& yoffset= screensize[1]-ysize	& END
6: BEGIN & xoffset = screensize[0]-xsize	& yoffset= yCenter-yHalfSize  	& END
7: BEGIN & xoffset = 0   					& yoffset= yCenter-yHalfSize  	& END
8: BEGIN & xoffset = xCenter-xHalfSize   	& yoffset= yCenter-yHalfSize  	& END
ELSE:
ENDCASE

WIDGET_CONTROL, id_widget, XOFFSET=xoffset, YOFFSET=yoffset

RETURN, 1
END

;**************************************************************************************************
;**************************************************************************************************
FUNCTION Interface_LabelInfo_ok, ev

	WIDGET_CONTROL, ev.top, GET_UVALUE=ptr_ok
	*ptr_ok = 'OK'
	WIDGET_CONTROL, ev.top, /DESTROY

END
;**************************************************************************************************
;**************************************************************************************************
FUNCTION Interface_LabelInfo_cancel, ev

	WIDGET_CONTROL, ev.top, GET_UVALUE=ptr_ok
	*ptr_ok = 'CANCEL'
	WIDGET_CONTROL, ev.top, /DESTROY

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION GetWidget_RelativePosition, id_widget, ID_PARENT=id_parent, POSITION=position, NOOUTSIDE=nooutside

; Relative position [x,y] in a widget, with respect to its parent

; Possible positions
;
;    _________
;  9|0   1   2|10
;   |         |
; 11|3   4   5|12
;   |         |
; 13|6   7   8|14
;    ---------
;

; nooutside : If the chosen position is outside the screen, returns, -1

IF N_ELEMENTS(position) NE 0 THEN pos = FIX(position) ELSE pos = 2


geom_parent = WIDGET_INFO(id_parent, /GEOMETRY)
geom_widget = WIDGET_INFO(id_widget, /GEOMETRY)
xsp = geom_parent.Scr_XSize
ysp = geom_parent.Scr_YSize
xpp = geom_parent.xoffset
ypp = geom_parent.yoffset

xsw = geom_widget.Scr_XSize
ysw = geom_widget.Scr_YSize

CASE pos OF
0 : BEGIN & xoff = xpp 					& yoff= ypp & END
9 : BEGIN & xoff = xpp - xsw			& yoff= ypp & END
1 : BEGIN & xoff = xpp - xsw/2 + xsp/2 	& yoff= ypp & END
10: BEGIN & xoff = xpp + xsp			& yoff= ypp & END
2 : BEGIN & xoff = xpp + xsp - xsw		& yoff= ypp & END

11: BEGIN & xoff = xpp - xsw			& yoff= ypp - ysw/2 + ysp/2 & END
3 : BEGIN & xoff = xpp 					& yoff= ypp - ysw/2 + ysp/2 & END
4 : BEGIN & xoff = xpp - xsw/2 + xsp/2 	& yoff= ypp - ysw/2 + ysp/2 & END
12: BEGIN & xoff = xpp + xsp			& yoff= ypp - ysw/2 + ysp/2 & END
5 : BEGIN & xoff = xpp + xsp - xsw 		& yoff= ypp - ysw/2 + ysp/2 & END

13: BEGIN & xoff = xpp - xsw			& yoff= ypp - ysw + ysp & END
6 : BEGIN & xoff = xpp 					& yoff= ypp - ysw + ysp & END

7:  BEGIN & xoff = xpp - xsw/2 + xsp/2 	& yoff= ypp - ysw + ysp & END
14: BEGIN & xoff = xpp + xsp			& yoff= ypp - ysw + ysp & END
8:  BEGIN & xoff = xpp + xsp - xsw 		& yoff= ypp - ysw + ysp & END

ELSE:
ENDCASE

opt_nooutside = KEYWORD_SET(nooutside)

IF KEYWORD_SET(nooutside) THEN BEGIN
	DEVICE, Get_Screen_Size = screenSize
	IF ((xoff+xsw) GT screenSize[0]) THEN RETURN,-1
	IF ((yoff+ysw) GT screenSize[1]) THEN RETURN,-1
	IF (xoff LT 0) THEN RETURN,-1
	IF (yoff LT 0) THEN RETURN,-1
ENDIF

RETURN, [xoff, yoff]

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION IdxPoints_ChangeYconvention, idx_array, DIMENSIONS=dimensions, ARRAY_XY=array_xy

; Change index coordinates convention from Y increasing (IDL convention) to y from up to down
; (YUP convention)

; index in IDL convention and IMAFEN convention = y*sizex + x

xsize= dimensions[0]
ysize= dimensions[1]

xarray     = idx_array MOD xsize
yarray_old = idx_array/xsize
yarray   = ysize-yarray_old-1
idx_new  = yarray*xsize + xarray;


xarray+=1 ; begin with one
yarray+=1
array_xy = [TRANSPOSE(xarray),TRANSPOSE(yarray)]

RETURN, idx_new

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION IdxPoints_FirstFromlowRes, MASK_IDX=mask_idx, DIMENSIONS=dimensions, ROI_RESOLUTION=roi_resolution,$
			MASK_OUTPUT=mask_output

; retorna un array de índices (x+y*xsize) de los primeros pixeles (esquina superior izquierda)
; en alta resolución, de los pixeles de baja resolución

idx_array = WHERE(mask_idx GE 0, ct)

xsize = dimensions[0]
ysize = dimensions[1]

IF roi_resolution EQ 1 THEN BEGIN
    mask_output = BYTARR(xsize, ysize)
	mask_output[idx_array]=1
	RETURN, idx_array
ENDIF

np_pixlr   = roi_resolution*roi_resolution
pos_corner = np_pixlr - roi_resolution ; esquina superior izquierda

max_val     = MAX(mask_idx)
idx_arr_sel = INTARR(max_val+1)

FOR i=0l, max_val DO BEGIN
	pos = WHERE(mask_idx EQ i, ct)
	IF ct NE np_pixlr THEN RETURN,-1
	idx_arr_sel[i] = pos[pos_corner]
ENDFOR

mask_output = BYTARR(xsize, ysize)
mask_output[idx_array]=1
mask_output[idx_arr_sel]=2 ; auxiliar result


RETURN, idx_arr_sel

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION  IdxPointsFromBox, POINT_INI=point_ini, SIZE_ROI=size_roi, DIMENSIONS=dimensions

; indices de los puntos (alta resolucion) a partir de un punto inicial y un box
; para hacer el caso de box cuadrado un caso particular del general (free)

mask = BYTARR(dimensions[0:1])

mask[point_ini[0]:point_ini[0]+size_roi[0]-1,point_ini[1]:point_ini[1]+size_roi[1]-1] = 1

points_idx = WHERE(mask EQ 1, ct)

RETURN, points_idx

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION IdxPointsLr_From_Indexes, IDX_POINTS=idx_points, ROI_RESOLUTION=roi_resolution, MIN_P=min_p,$
	MAX_P=max_p, DIMENSIONS=dimensions, MASK_LR=mask_lr, MASK_HR=mask_hr, ERROR_ID=error_id, ERROR_STR=error_str

; A partir de los puntos (alta resolución) la resolución de la ROI devuelve los indices en baja resolucion,
; el minimo y maximo (alta resolución) y la máscara, (baja resolución)

; inputs:

; IDX_POINTS     =puntos alta resolución
; DIMENSIONS     =dimensiones
; ROI_RESOLUTION= resolución de los vóxeles

	error_id = 0
	xsize = dimensions[0]
	ysize = dimensions[1]

	mask_hr = BYTARR(xsize, ysize)
	mask_hr[idx_points]=1

	limits = IniPointfromIndexes(IDX_POINTS=idx_points, DIMENSIONS=dimensions, ROI_RESOLUTION=roi_resolution,$
		MIN_P=min_p, MAX_P=max_p)

	mask_hr = mask_hr[min_p[0]:max_p[0],min_p[1]:max_p[1]]

	xsize_subim = (max_p[0]-min_p[0]+1)
	ysize_subim = (max_p[1]-min_p[1]+1)

	IF  xsize_subim MOD roi_resolution NE 0 THEN BEGIN
		error_id = 1 &	RETURN, -1
	ENDIF
	IF  ysize_subim MOD roi_resolution NE 0 THEN BEGIN
		error_id = 1 &	RETURN, -1
	ENDIF

	xsize_subim_lr = xsize_subim/roi_resolution
	ysize_subim_lr = ysize_subim/roi_resolution

	IF roi_resolution GT 1 THEN BEGIN
		mask_lr  = REBIN(FLOAT(mask_hr), xsize_subim_lr, ysize_subim_lr) GE 0.5
		IF TOTAL(mask_lr,/INTEGER) EQ 0 THEN BEGIN ; avoid total equal cero
			mask_lr = REBIN(FLOAT(DILATE(mask_hr, get_disc(4))), xsize_subim_lr, ysize_subim_lr) GE 0.5
		ENDIF
	ENDIF ELSE mask_lr = mask_hr

	;IF roi_resolution GT 1 THEN BEGIN
	;	mask_lr  = REBIN(mask_hr, xsize_subim_lr, ysize_subim_lr)
	;ENDIF ELSE mask_lr = mask_hr

	pos_idx_sb = WHERE(mask_lr NE 0, ct)
	IF ct EQ 0 THEN BEGIN
		error_id = 2 & error_str = 'Small ROI, no voxels inside' & RETURN, -1
	ENDIF

	RETURN, pos_idx_sb

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Imag2ThreeChannels, image, RED_VALUES=red_Values, GREEN_VALUES=green_values, BLUE_VALUES=blue_values

IF SIZE(image,/TYPE) NE 1 THEN RETURN, -1;byte value

n_dim = SIZE(image, /N_DIMENSIONS)
CASE n_dim OF
1 : BEGIN
	xsize = (SIZE(image, /DIMENSIONS))[0]
	ysize = 1L
END
2 : BEGIN
	xsize = (SIZE(image, /DIMENSIONS))[0]
	ysize = (SIZE(image, /DIMENSIONS))[1]

END
ELSE: RETURN, -1
ENDCASE

imag_color = BYTARR(xsize,ysize,3)

imag_color[*,*,0] = red_Values[image]
imag_color[*,*,1] = green_values[image]
imag_color[*,*,2] = blue_values[image]

RETURN, imag_color

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION ImagFromArray, array, DIMENSIONS=dimensions, $
	IDX_POINTS = idx_points, PINI=pini, MASK=mask, IDX_NOVALID=idx_novalid

type = SIZE(array, /TYPE)
xsize = dimensions[0]
ysize = dimensions[1]

n_points = N_ELEMENTS(array)
IF n_points NE N_ELEMENTS(idx_points) THEN RETURN, -1

idxarr_x = idx_points MOD xsize
idxarr_y = idx_points / xsize

image = MAKE_ARRAY(xsize, ysize, TYPE=type)
mask  = BYTARR(xsize,ysize)
mask[idx_points]  = 1
image[idx_points] = REFORM(array)

IF N_ELEMENTS(idx_novalid) NE 0 THEN BEGIN
	mask[idx_points[idx_novalid]] = 0
ENDIF

RETURN, image

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION ImIdxs_From_Indexes, MASK_LR=mask_lr, DIMENSIONS=dimensions, ROI_RESOLUTION=roi_resolution, $
	MIN_P = min_p, IDX_POINTS=idx_points

xsize = dimensions[0]
ysize = dimensions[1]

xsize_roi = (SIZE(mask_lr, /DIMENSIONS))[0]
ysize_roi = (SIZE(mask_lr, /DIMENSIONS))[1]

points_idx = WHERE(mask_lr NE 0, ct)
IF ct EQ 0 THEN RETURN, -1

mask_idx = INTARR(xsize_roi, ysize_roi)-1 ; valores -1 no tienen point_idx
mask_idx[points_idx] = INDGEN(N_ELEMENTS(points_idx))


IF roi_resolution GT 1 THEN BEGIN
	mask_idx = REBIN(mask_idx, xsize_roi*roi_resolution, ysize_roi*roi_resolution, /SAMPLE)
	; sample mask_idx
ENDIF

mask_tot = INTARR(xsize,ysize)-1
max_p = min_p + [xsize_roi,ysize_roi]*roi_resolution - 1
IF max_p[0] GE xsize THEN RETURN,-1
IF max_p[1] GE ysize THEN RETURN,-1
mask_tot[min_p[0]:max_p[0] , min_p[1]:max_p[1]] = mask_idx

idx_points = WHERE(mask_tot NE -1, ct)

RETURN ,mask_tot

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION IniPointfromIndexes, IDX_POINTS=idx_points, DIMENSIONS=dimensions, MIN_P=min_p, MAX_P=max_p,$
ROI_RESOLUTION=roi_resolution

xsize = dimensions[0]
ysize = dimensions[1]

IF N_ELEMENTS(roi_resolution) EQ 0 THEN roi_resolution=1l

idxarr_x = idx_points MOD xsize
idxarr_y = idx_points / xsize

min_x = MIN(idxarr_x, MAX=max_x)
min_y = MIN(idxarr_y, MAX=max_y)

roi_size = [max_x-min_x+1,max_y-min_y+1]
roi_size_ini = roi_size

IF roi_size[0] LT roi_resolution THEN roi_size[0] = roi_resolution
IF roi_size[1] LT roi_resolution THEN roi_size[1] = roi_resolution

IF roi_size[0] MOD roi_resolution NE 0 THEN roi_size[0] = roi_size[0] + (roi_resolution-(roi_size[0] MOD roi_resolution))
IF roi_size[1] MOD roi_resolution NE 0 THEN roi_size[1] = roi_size[1] + (roi_resolution-(roi_size[1] MOD roi_resolution))

IF roi_size[0] GT xsize THEN roi_size[0]-=roi_resolution
IF roi_size[1] GT ysize THEN roi_size[1]-=roi_resolution

diff = roi_size-roi_size_ini

IF diff[0] NE 0 THEN BEGIN ; distribuye la desigualdad equitativamente
	min_x-=diff[0]/2
	max_x-=diff[0]/2
	IF (max_x-min_x+1) NE roi_size[0] THEN max_x+=(roi_size[0]-(max_x-min_x+1))
ENDIF
IF diff[1] NE 0 THEN BEGIN ; distribuye la desigualdad equitativamente
	min_y-=diff[1]/2
	max_y-=diff[1]/2
	IF (max_y-min_y+1) NE roi_size[1] THEN max_y+=(roi_size[1]-(max_y-min_y+1))
ENDIF

;------------------------------------------
; Cuidado con los bordes
IF max_x GT (xsize-1) THEN BEGIN
	max_x-=max_x-(xsize-1)
	min_x-=max_x-(xsize-1)
ENDIF
IF max_y GT (ysize-1) THEN BEGIN
	max_y-=max_y-(ysize-1)
	min_y-=max_y-(ysize-1)
ENDIF
IF min_x LT 0  THEN BEGIN
	max_x-=min_x
	min_x=0
ENDIF
IF min_y LT 0 THEN BEGIN
	max_y-=min_y
	min_y =0
ENDIF
;-------------------------------------------
size_roi = [max_x-min_x+1,max_y-min_y+1]

; solo puede ser multiplo
IF (max_x-min_x+1) MOD roi_resolution NE 0 THEN RETURN, -1
IF (max_y-min_y+1) MOD roi_resolution NE 0 THEN RETURN, -1

; Ahora evita ROIs de solo un punto en X e Y (en baja resolucion)
IF (max_x - min_x+1) EQ roi_resolution THEN BEGIN
	IF max_x LT (xsize-roi_resolution) THEN max_x+=roi_resolution ELSE min_x-=roi_resolution
ENDIF
IF (max_y - min_y+1) EQ roi_resolution THEN BEGIN
	IF max_y LT (ysize-roi_resolution) THEN max_y+=roi_resolution ELSE min_y-=roi_resolution
ENDIF

min_p = [min_x, min_y]
max_p = [max_x, max_y]

RETURN, [[min_x, min_y],[max_x, max_y]]

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
PRO Interface_datain_destroy, ev

	WIDGET_CONTROL, ev, GET_UVALUE=id
	WIDGET_CONTROL,  id.wbase_main,  /DESTROY

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<++>

PRO Interface_datain_generic_event, ev

WIDGET_CONTROL, ev.top, GET_UVALUE=id
WIDGET_CONTROL, ev.id,  GET_UVALUE=uval

opt_preview = 1

FOR i=0l, N_ELEMENTS(id.wtext_id)-1 DO BEGIN
	WIDGET_CONTROL, id.wtext_id[i], GET_VALUE=text_Arr
	*(id.ptr_text[i]) = text_arr[0]
	PRINT, text_arr
ENDFOR
WIDGET_CONTROL, id.wbase_main, SET_UVALUE= id

IF N_ELEMENTS(uval) GT 0 THEN BEGIN
CASE uval OF
'EXIT1':BEGIN
	Interface_datain_generic_destroy, ev.top
END
'EXIT2':BEGIN
	FOR i=0l, N_ELEMENTS(id.wtext_id)-1 DO BEGIN
		WIDGET_CONTROL, id.wtext_id[i], SET_VALUE=''
		*(id.ptr_text[i]) = ''
	ENDFOR
	WIDGET_CONTROL, id.wbase_main, SET_UVALUE= id
	Interface_datain_generic_destroy, ev.top

END
ELSE:
ENDCASE
ENDIF

END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_datain_generic_destroy, ev

	WIDGET_CONTROL, ev, GET_UVALUE=id
	WIDGET_CONTROL,  id.wbase_main,  /DESTROY

END
;**************************************************************************************************
;**************************************************************************************************


FUNCTION Interface_datain_generic, GROUP_LEADER=group_leader, XSIZE=xsize, YSIZE=ysize, POSITION=position, $
    TITLE=title, ARRAY_QUESTIONS=array_questions, RESULT=result, INITIAL_RESULT=initial_result, ROW=row

IF KEYWORD_SET(row) THEN opt_row=1 ELSE opt_row=0

tlb_frame_attr=1

n_questions = N_ELEMENTS(array_questions)
IF N_ELEMENTS(initial_result) NE n_questions THEN ini_result=STRARR(n_questions) ELSE $
	ini_result = STRING(initial_result)

IF N_ELEMENTS(group_leader) NE 0 THEN BEGIN
    base_floating = WIDGET_BASE(GROUP_LEADER=group_leader, FLOATING=1, UVALUE='INFO_FLOATING', $
       TITLE=title, TLB_FRAME_ATTR=tlb_frame_attr, MAP=opt_map, COLUMN=1)
ENDIF ELSE BEGIN
    base_floating = WIDGET_BASE(UVALUE='INFO_FLOATING', $
       TITLE=title, TLB_FRAME_ATTR=tlb_frame_attr, MAP=opt_map, COLUMN=1)
ENDELSE

base_0 = WIDGET_BASE(base_floating, COLUMN=opt_row EQ 0, ROW=opt_row, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
arr_bases = LONARR(n_questions)
arr_textid  = LONARR(n_questions)
FOR i=0, n_questions-1 DO BEGIN
	arr_bases[i]  = WIDGET_BASE(base_0, ROW=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
	label         = WIDGET_LABEL(arr_bases[i], VALUE=array_questions[i])
	arr_textid[i] = WIDGET_TEXT(arr_bases[i], /EDITABLE, VALUE=ini_result[i], XSIZE=8)
ENDFOR
base_1 = WIDGET_BASE(base_floating, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
	bttn_close  = WIDGET_BUTTON(base_1, VALUE=' OK     ', UVALUE='EXIT1')
	bttn_cancel = WIDGET_BUTTON(base_1, VALUE=' Cancel ', UVALUE='EXIT2')


IF N_ELEMENTS(position) EQ 1 AND N_ELEMENTS(group_leader) NE 0 THEN BEGIN
	pos = GetWidget_RelativePosition(base_floating, ID_PARENT=group_leader, POSITION=position)
	WIDGET_CONTROL, base_floating, XOFFSET=pos[0], YOFFSET=pos[1]
ENDIF ELSE BEGIN
	 pos = GetWidget_PositionInScreen(base_floating,CENTER=1)
ENDELSE

st = { $
	wbase_main : base_floating,$
	wbttn_close : bttn_close,$
	wtext_id   : arr_textid,$
	ptr_text   : PTRARR(n_questions, /ALLOCATE_HEAP),$
	last : 0l }

FOR i=0,n_questions-1 DO BEGIN
	*(st.ptr_text[i]) = ''
ENDFOR

WIDGET_CONTROL, base_floating, SET_UVALUE=st
WIDGET_CONTROL, base_floating, /REALIZE

XMANAGER,  'Interface_datain_generic',  base_floating, NO_BLOCK=1, MODAL=1, $
    EVENT_HANDLER='Interface_datain_generic_event', CLEANUP='Interface_datain_generic_Destroy'

PRINT, 'OK'

result = STRARR(n_questions)

FOR i=0l, n_questions-1 DO BEGIN
	IF N_ELEMENTS(*(st.ptr_text[i])) EQ 1 THEN $
		result[i] = *(st.ptr_text[i]) $
	ELSE result[i]=''
ENDFOR

PTR_FREE, st.ptr_text
; Interfaz modal de entrada de datos

RETURN, result

END

;<-->
;<++>

FUNCTION ParametricCp_to_CurveCp ,MODEL=model, ST_AIF=st_aif, PARAMS=params, ARR_TIME=arr_time, INJECTION_FRAME=injection_frame, $
			HAEMATOCRIT=haematocrit, ERROR_ID=error_id, DOSE=dose

error_id=0

arr_time_inj = DOUBLE(arr_time[injection_frame:N_ELEMENTS(arr_time)-1]) - arr_time[injection_frame]  ; origin of times in the injection frame

;------------------------------------------------
IF N_ELEMENTS(st_aif) NE 0 THEN BEGIN
	IF SIZE(st_aif,/TYPE) EQ 8 THEN BEGIN
		option_st_aif = 1
		str_model = STRUPCASE(st_aif.model)
	ENDIF ELSE BEGIN
		error_id=1 & RETURN, -1
	ENDELSE
ENDIF ELSE BEGIN
	IF N_ELEMENTS(model) NE 0 THEN BEGIN
		str_model = STRUPCASE(model)
		option_st_aif = 0
	ENDIF ELSE BEGIN
		error_id=2 & RETURN, -1
	ENDELSE
ENDELSE
;-----------------------------------------------
arr_models = ['BIEXPONENTIAL', 'ORTON','SCHABEL-SIMPLIFIED', 'SCHABEL-STANDARD', 'PARKER']

num_model = (WHERE(str_model EQ arr_models,ct))[0]

CASE num_model OF
	0 : BEGIN  ;'BIEXPONENTIAL':

		arr_cp = Model_Cp_Biexponential(DOSE=dose, ARR_TIME=arr_time_inj, PARAMS=params, STRUCT=st_aif)
		arr_cp = [0, arr_cp[0:N_ELEMENTS(arr_cp)-2]]; retrasa para compatibilidad
	END
	1 : BEGIN  ;'ORTON':
		arr_cp = Model_Cp_Orton2(ARR_TIME_MIN=arr_time_inj, PARAMS=params, STRUCT=st_aif)
	END
	2:  BEGIN ;'SCHABEL_SIMPLIFIED':
		arr_cp = Model_Cp_SchabelSimplified(ARR_TIME_MIN=arr_time_inj, PARAMS=params, STRUCT=st_aif)
	END
	3:  BEGIN ;'SCHABEL':
		;params  = [ A[1]/A[0], A[2]/A[0], A[3]/A[0],  alpha[0],$
;		    tau[0], tau[1], delay_t[0], delta_t[0], delta_t[1]-delta_t[0], delta_t[2]-delta_t[0]]
		params_10_in = params[1:10]
		params_10_in[0:2]/=params[0]
		params_10_in[8:9]-=params[8]
		arr_cp = Model_Cp_Schabel(ARR_TIME_MIN=arr_time_inj, DELTA_T_MIN=delta_t_min, $
			PARAMS_10_IN =params_10_in, A0=params[0])
	END
	4 : BEGIN; 'PARKER':
		arr_cp = Model_Cp_Parker(ARR_TIME_MIN=arr_time_inj, PARAMS=params, STRUCT=st_aif)
	END
	ELSE : BEGIN
		error_id = 5
		RETURN, -1
	END
ENDCASE

IF injection_frame GT 0 THEN BEGIN
	arr_zeros = REPLICATE(0, injection_frame)
	RETURN, [arr_Zeros, arr_cp]
ENDIF ELSE BEGIN
	RETURN, arr_cp
ENDELSE

END

;**************************************************************************************************
;**************************************************************************************************

PRO DCEMRI_WelcomeWidget, OK=ok

DEVICE, DECOMPOSED=0
file_logo1 = '.\Icons\Logobit_ciber.png'
file_logo2 = 'Logobit_ciber.png'
file_logo3 = '..\Icons\Logobit_ciber.png'
IF FILE_TEST(file_logo1) THEN image=READ_PNG(file_logo1, R,G,B)
IF FILE_TEST(file_logo2) THEN image=READ_PNG(file_logo2, R,G,B)
IF FILE_TEST(file_logo3) THEN image=READ_PNG(file_logo3, R,G,B)


IF N_ELEMENTS(image) LE 1 THEN BEGIN
	image = BYTARR(3,2,2)+255
	wait_temp = 0
ENDIF ELSE wait_temp = 1

DEVICE, Get_Screen_Size=screenSize
fullScrXsize=screensize[0]
fullScrYsize=screensize[1]
s = SIZE(Image, /DIMENSIONS)
nx = s[1]
ny = s[2]

splashXoffset = (fullScrXsize-nx)/2
splashYoffset = (fullScrYsize-ny)/2

portada_Base = WIDGET_BASE( tlb_frame_attr=31,$
                           xoffset=splashXoffset, yoffset=splashYoffset,xsize=nx+3,ysize=ny+3)
splashDraw = WIDGET_DRAW( portada_Base, xsize=nx, ysize=ny, /frame)

WIDGET_CONTROL, portada_Base, MAP=0
WIDGET_CONTROL, portada_Base, /REALIZE

WIDGET_CONTROL, splashDraw, GET_VALUE=splashDrawID
WSET, splashDrawID
TV, Image,TRUE=1
Image = 0                 ;Don't need it anymore
WIDGET_CONTROL, portada_Base, MAP=1

message_text = Get_license_array()

WAIT, wait_temp

	result = Interface_LabelInfo(GROUP_LEADER=portada_Base, XSIZE=180, YSIZE=N_ELEMENTS(message_text), POSITION=position,$
    	TITLE='Disclaimer', LABEL=Message_Text, NO_CLOSE=no_close, MODAL=1, FONT_SIZE=font_size, FONT_TYPE=font_type, FONT_EFFECT=font_effect,$
   		HORIZONTAL=horironzal, OKCANCEL=1)

	if result eq 'OK' then OK=1 ELSE OK=0


WIDGET_CONTROL, portada_Base, /DESTROY
DEVICE, DECOMPOSED=1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION DCEMRI_ReadConfig, ARRAY_IN=array_in, ARRAY_OUT=array_out, $
	ERROR_ID=error_id, ERROR_STR=error_str, ST_AUX=st_aux, ST_AIF=st_aif

error_id   =  0
error_str  = ''

n_lines_min = 10

n_lines = N_ELEMENTS(array_in)
strarr_lines = array_in;
IF n_lines LT n_lines_min THEN BEGIN & error_id = 4 & error_str = 'Error #3, not valid config file' & RETURN, -1 & ENDIF

;------------------------------------------------------
strarr_values = ReadStrarr_obtainvalues(strarr_lines[0:n_lines-1], STR_SEPARATOR='=', ERROR_ID=err2, ERROR_STR=error_str)
IF err2 NE 0 THEN BEGIN & error_id =7l &  RETURN, -1 & ENDIF
strarr_values[0,*] = STRLOWCASE(strarr_values[0,*])
;------------------------------------------------------

strarr_data = [$
	'path_data',  $
	'path_result',$
	'path_t1map', $
	'path_aif',   $

	'size_window',    $
	'size_roi',$
	'type_roi',$
	'resolution_roi', $
	'interpolation',  $


	't10_tissue',$
	't10_blood', $
	'repetition_time',$
	'flip_angle',     $
	'frame_period',   $

	'frame_injection',$
	'injected_dose',  $
	'relaxivity',  $
	'haematocrit', $


	'size_im',$
	'nslices',$
	'nframes',$
	'nframes_iauc',$
	;'nframes_stats',$

	'aif_m1',$
	'aif_m2',$
	'aif_a1',$
	'aif_a2',$

	'ktrans_ini',$
	'kep_ini',$
	'vp_ini' ,$

	'ah_ini'  ,$
	;'bh_ini'   ,$
	'keph_ini' ,$
	'kelh_ini' ,$
	'tau_ini'  ,$

	'kepl_ini',$
	'sl_ini' ,$

	'ktrans_rr' ,$
	've_rr', $

	;'data_err',$
	'tau_infusion' ]

strarr_Data = STRLOWCASE(strarr_data)
;---------------------------------------------
n_data  = N_ELEMENTS(strarr_Data)
err=0
ascii_min = BYTE('0')
ascii_max = BYTE('9')

n_campos=0l
FOR i=0l, n_data-1 DO BEGIN
	pos = WHERE(STRMATCH(strarr_values[0,*], strarr_data[i]) EQ 1, ct)
	IF ct NE 1 THEN BEGIN
		err = -1
		error_str = 'Keyword not found: ' + strarr_data[i]
		BREAK
	ENDIF
	var_temp = (STRTRIM((strarr_values[1,pos])[0],2))[0]

	 CASE strarr_data[i] OF
	 	'path_data'     : path_data   = var_temp[0]
	 	'path_result'   : path_result = var_temp[0]
	 	'path_t1map'    : path_t1map  = var_temp[0]
	 	'path_aif'		: path_aif    = var_temp[0]

	 	'size_window'   : size_window = FIX(var_temp)
	 	'size_roi'      : size_roi = ReadArray_numbers(var_temp, /INTEGER)
	 	'resolution_roi': resolution_roi = FIX(var_temp)
	 	'interpolation' : interpolation = (STRTRIM(STRUPCASE(var_temp),2) EQ 'YES')

		'injected_dose'   : injected_dose = FLOAT(var_temp)
		'frame_injection' : frame_injection_begin1 = FIX(var_temp)
		'frame_period'    : frame_period_seconds = FLOAT(var_temp)
		't10_tissue'      : T10_tissue = (FLOAT(var_temp))[0]
		't10_blood'       : T10_blood = (FLOAT(var_temp))[0]
		'repetition_time' : TR = (FLOAT(var_temp))[0]
		'relaxivity'      : relaxivity = (FLOAT(var_temp))[0]
		'haematocrit'     : haematocrit = (FLOAT(var_temp))[0]

		'flip_angle'	  : flip_angle_degrees = (FLOAT(var_temp))[0]

		'size_im'		  : size_im = ReadArray_numbers(var_temp, /INTEGER)
		'nslices'         : nslices = FIX(var_temp)
		'nframes'         : nframes = FIX(var_temp)
		'nframes_iauc'    : frames_for_iauc = FIX(var_temp)
		;'nframes_stats'   : nframes_stats = ReadArray_numbers(var_temp, /INTEGER)

		'aif_m1' : m1 = FLOAT(var_temp)
		'aif_m2' : m2 = FLOAT(var_temp)
		'aif_a1' : a1 = FLOAT(var_temp)
		'aif_a2' : a2 = FLOAT(var_temp)

		'ktrans_ini' : ktrans_ini = FLOAT(var_temp)
		'kep_ini'    : Kep_ini    = FLOAT(var_temp)
		'vp_ini'     : vp_ini     = FLOAT(var_temp)

		'ah_ini'   : Ah_ini   = FLOAT(var_temp)
		;'bh_ini'   : Bh_ini   = FLOAT(var_temp)
		'keph_ini' : keph_ini = FLOAT(var_temp)
		'kelh_ini' : kelh_ini = FLOAT(var_temp)
		'tau_ini'  : tau_ini  = FLOAT(var_temp)

		'kepl_ini' : kepl_ini =  FLOAT(var_temp)
		'sl_ini'   : sl_ini =  FLOAT(var_temp)

		'tau_infusion' : tau_infusion = FLOAT(var_temp)
		;'data_err' : data_err = FLOAT(var_temp)

		'ktrans_rr' : ktrans_rr = FLOAT(var_temp)
		've_rr'     : ve_rr =  FLOAT(var_temp)

		'type_roi' : BEGIN
				str_type_roi = STRUPCASE(var_temp)
				CASE str_type_roi OF  ;[0, cuadrada, 1: toda la imagen, 2: libre]
				'BOX': type_roi = 0l
				'FULL':type_roi = 1l
				'FREE':type_roi = 2l
				ELSE : BEGIN
					error_id = 18l
					error_str = 'not valid type of ROI : ' + str_type_roi
					RETURN,-1
				END
				ENDCASE
		END
	ELSE : BEGIN
		error_id = 8l
		error_str = 'not valid keyword : ' + strarr_data[i]
		RETURN,-1
	END
	ENDCASE
    n_campos++
ENDFOR
IF n_campos NE n_data THEN BEGIN & error_id = 9l & RETURN,-1 & ENDIF

;IF N_ELEMENTS(nframes_stats) NE 2 THEN BEGIN
;	error_id = 10l & error_str = 'NFRAMES_STATS requires two parameters in format [a,b]' & RETURN,-1
;ENDIF
IF N_ELEMENTS(size_roi) NE 2 THEN BEGIN
	error_id = 11l & error_str = 'SIZE_ROI requires two parameters in format [a,b]' & RETURN,-1
ENDIF
IF N_ELEMENTS(size_im) NE 2 THEN BEGIN
	error_id = 12l & error_str = 'SIZE_IM requires two parameters in format [a,b]' & RETURN,-1
ENDIF

str_comma =','

st_param = {$
;---------------------------------------------
	dose    		 : injected_dose,   $
	nframe_injection : frame_injection_Begin1-1, $
	frame_period     : frame_period_seconds/60d, $
	T10_tissue : T10_tissue, $  ; t10 por defecto cuando no se carga imagen, (en segundos)
	T10_blood  : T10_blood, $  ; t10 por defecto cuando no se carga imagen, (en segundos)
	TR  : TR , $
	flip_angle_degrees  : flip_angle_Degrees, $
	haematocrit         : haematocrit,$
	r1     		     	: relaxivity, $
	nframes             : nframes,    $
	nframes_auc         : frames_for_iauc, $
	nframes_stats       : [0,0l],          $
	ktrans_ini : ktrans_ini, $
	kep_ini    : Kep_ini ,   $
	vp_ini     : vp_ini , $
	Ah_ini   : Ah_ini   , $
	Bh_ini   : 0.0   ,    $
	keph_ini : keph_ini , $
	kelh_ini : kelh_ini,  $
	tau_ini  : tau_ini  , $
	kepl_ini : kepl_ini , $
	sl_ini   : sl_ini  ,  $
	tau_infusion : tau_infusion, $
	;data_err : data_err,$
	nframes_adjust : 0l, $ $
	ktrans_rr : ktrans_rr, $
	ve_rr     : ve_rr,      $
	kep_rr    : ktrans_rr/ve_rr    $
}
;---------------------------------------------
st_param.nframes_adjust = st_param.nframes-st_param.nframe_injection;
st_param.nframes_stats  = [0, st_param.nframes_adjust-1]
;---------------------------------------------
st_aux = { $
	sizeview : size_window,$
	size_im  : size_im,$
	nframes  : st_param.nframes,$
	nslices  : nslices,$
	size_roi : size_roi,$
	type_roi : type_roi,$
	np_lr    : resolution_roi,  $
	path_data   : path_data,    $
	path_result : path_result , $
	path_t1map  : path_t1map,   $
	path_aif    : path_aif,     $
	flag_interpolation : interpolation}
;---------------------------------------------
st_aif = { $
	model : 'biexponential',$
	m1 : m1,  $
	m2 : m2,  $
	a1 : a1,  $
	a2 : a2  }
;---------------------------------------------
error_str = ''
array_out = REFORM(strarr_lines, 1, N_ELEMENTS(strarr_lines))

RETURN, st_param

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION DCEMRI_InitialParameters, ST_AIF=st_aif

	nframes             = 50
	nframe_injection    = 2l  ;Cuidado, empieza en cero (aqui un valor menor que el que sale en pantalla)
	;nframe_injection   = 10l  ;Cuidado, empieza en cero (aqui un valor menor que el que sale en pantalla)
	nframes_adjust= nframes-nframe_injection
	nframes_auc   = 3l   ; frames sobre los que se mide la IAUC
	nframes_stats = [0,nframes_adjust-1]; Rango de frames sobre los que se mide la estadistica de RMSE y ME

	st_param = {$

		; Tofts parameters
		ktrans_ini: 0.03,$ ; min-1
		Kep_ini   : 0.18,$ ; min-1
		vp_ini    : 0.0, $ ; volume fraction
		; hoffmann parameters
		Ah_ini    :   1d, $
		;Bh_ini   :   0d, $;
		keph_ini  : 0.05, $
		kelh_ini  : 0.001,$
		tau_ini   : 0d,   $
		; Larsson parameters
		kepl_ini : 0.05,  $
		sl_ini   : 1d,    $

		tau_infusion : 0d,   $ ; infusion time (if perfect bolus injection: tau infusion=0)
		tr   		: 0.200, $ ; repetition time (seconds)
		t10_tissue  : 1.0,   $ ; relaxation time (seconds)
		t10_blood   : 1.0,   $ ; relaxation time in blood (seconds)
		dose 		: 0.2,$    ; mmol/Kg animal ;          (Dato de correo)
		r1 			: 3.6,$    ; (mM s-1)                  (Dato de correo)
		flip_angle_Degrees : 90.0, $
		haematocrit : 0.45, $


		frame_period     : 40.0/60,$ ; 25.6/60,  $ ; en minutos
		nframes          : nframes, $
		nframe_injection : nframe_injection,$
		nframes_adjust   : nframes_adjust, $
		nframes_auc      : nframes_auc, $
		nframes_stats    : nframes_stats[0:1],$

		ktrans_rr : 0.20,$
		ve_rr     : 0.20,$
		kep_rr    : 0.0  $

		;data_err : 6.0/100 $ ; Std dev of noise in acquisitions
	}
	st_param.kep_rr = st_param.ktrans_rr/st_param.ve_rr

	;--------------------------------------------------------------

	m1 = 2.810 ; en minutos  (0.8 min-1)  ; 2.70 min-1
	m2 = 0.023 ; en minutos  (0.021 min-1)
	a2 = (1/0.043)/(1.5833+1) ; kg/litro
	a1 = a2*1.5833            ; kg/litro

	st_aif = {$
		; aif parameters (Tofts, larsson)
		model: 'biexponential',$
		m1 : m1,$
		m2 : m2,$
	   	a2 : a2,$
		a1 : a1	}
	;--------------------------------------------------------------

	RETURN, st_param

END
;**************************************************************************************************
;**************************************************************************************************

FUNCTION DCEMRI_ConfigModelParameters, modelnumber, STR_MODEL=str_model, ERROR_ID=error_id, ARR_MODELS=arr_models, STRARR_MODELS=strarr_models

error_id    = 0

arr_models    = ['CURVE', 'HOFFMANN', 'TOFTS', 'TOFTS.EXT', 'LARSSON', 'RR']
strarr_models = ['Curve parameters','Hoffmann model','Tofts model','Extended Tofts model', 'Larsson model','RR model']
arr_nparam    = [3,6,5,6,4,5] ; Number of parameters by model

arr_typemodel  = [1,2,3,3,4,5] ; Tofts y extended totfs share the same type model

arr_using_T10  = [0,0,1,1,0,1] ; The model needs T1 map
arr_using_Ct   = [0,0,1,1,0,1] ; the model needs the tissue concentration (not only signal intensity)

IF N_ELEMENTS(modelnumber) EQ 0 THEN BEGIN
	nmodel = WHERE(str_model EQ arr_models, npos)
	IF npos NE 1 THEN RETURN,-1
ENDIF ELSE BEGIN
	nmodel = modelnumber[0]
ENDELSE

model      = arr_models[nmodel]
str_model  = strarr_models[nmodel]
nparam     = arr_nparam[nmodel]
typemodel  = arr_typemodel[nmodel]
using_t10  = arr_using_t10[nmodel]
using_Ct   = arr_using_ct[nmodel]

strarr_idsinf = STRARR(4, nparam)
nl = 0l

CASE model OF
'CURVE' : BEGIN

	; CHAPUZA MINDT
	opt_mindt = 0
	IF opt_mindt EQ 0 THEN BEGIN

		strarr_droplparams = ['RCE','IAUC', 'TTM'] ;These are the non parametric data
		strarr_ids = ['RCE','IAUC', 'TTM']
		sd_valid = [0,0,0]
		strarr_idsinf[0:2,*] = [  $
			['RCE =',     ' %',      'RCE (%)'    ],$
			['IAUC=',     ' a.u',    'IAUC (a.u)' ],$
			['TTM =',     ' min',    'TTM (min)'  ]]
		strarr_idsinf[3,nl]   = 'Dynamic constrast enhancement'
		strarr_idsinf[3,nl++] = 'Initial Area under curve'
		strarr_idsinf[3,nl++] = 'Time to Max'
	ENDIF ELSE BEGIN

		strarr_droplparams = ['rCBF','rCBV', 'rMTT'] ;These are the non parametric data
		strarr_ids = ['RCBF','RCBV', 'RMTT']
		sd_valid = [0,0,0]
		strarr_idsinf[0:2,*] = [  $
			['rCBF=',    ' a.u',   'rCBF (a.u)' ],$
			['rCBV=',    ' a.u',   'rCBV (a.u)' ],$
			['rMTT=',   ' a.u',    'rMTT (a.u)' ]]
		strarr_idsinf[3,nl]   = 'Relative cerebral blood flood'
		strarr_idsinf[3,nl++] = 'Relative cerebral blood volume'
		strarr_idsinf[3,nl++] = 'Relative mean transit time'
	ENDELSE

END
'HOFFMANN'   : BEGIN

	strarr_droplparams = ['A.kep', 'Ah', 'kep', 'kel', 'RMSE','ME'] ; parameters of Hoffmann-Brix models
	strarr_ids = ['AKEP', 'AH', 'KEP', 'KEL', 'RMSE','ME']
	sd_valid = [1,1,1,1,0,0]
	strarr_idsinf[0:2,*] = [ $
		['A.kep(h) =',' min-1',  'A*kep (min-1)'],$
		['A(h) =',      ' n.d',  'A (n.d)'    ],$
		['kep(h) =',  ' min-1',  'kep (min-1)'],$
		['kel(h) =',  ' min-1',  'kel (min-1)'],$
		['RMSE =',  ' RCE (%)',  'RMSE (%)'],$
		['ME =',    ' RCE (%)',  'ME (%)']]
	strarr_idsinf[3,nl++] = 'Amplitude * kep (Hoffmann model)' ; A*kep
	strarr_idsinf[3,nl++] = 'Amplitude (Hoffmann model)' ; A
	strarr_idsinf[3,nl++] = 'Rate wash-in constant between EES and blood plasma (Hoffmann model)' ; kep (h)
	strarr_idsinf[3,nl++] = 'Rate wash-out constant between EES and blood plasma (Hoffmann model)'
	strarr_idsinf[3,nl++] = 'Root mean squared error (RMSE) of Hoffmann fitting'
	strarr_idsinf[3,nl++] = 'Mean absolute error (ME) of Hoffmann fitting'


END
'TOFTS'     : BEGIN
	strarr_droplparams = ['Ktrans', 'kep', 've', 'RMSE','ME']       ; parameters of Tofts models
	strarr_ids = ['KTRANS', 'KEP', 'VE', 'RMSE','ME']
	sd_valid = [1,1,1,0,0]
	strarr_idsinf[0:2,*] = [ $
		['Ktrans =',  ' min-1',  'Ktrans (min-1)'],$
		['kep =',     ' min-1',  'kep (min-1)'   ],$
		['ve =',      ' ',       've'            ],$
		['RMSE =',  ' mM (Ct)',  'RMSE  (mM)'  ],$
		['ME =',    ' mM (Ct)',  'ME (mM)']]
	strarr_idsinf[3,nl++] = 'Volume transfer constant between blood plasma and EES (Tofts model)' ; KTRANS
	strarr_idsinf[3,nl++] = 'Rate constant between EES and blood plasma (Tofts model)'
	strarr_idsinf[3,nl++] = 'Extracellular-extravascular instertitial space (EES) per unit volume of tissue (Tofts model)'
	strarr_idsinf[3,nl++] = 'Root mean squared error (RMSE) Tofts fitting'
	strarr_idsinf[3,nl++] = 'Mean absolute error (ME) of Tofts fitting'

END
'TOFTS.EXT' : BEGIN
	strarr_droplparams = ['Ktrans', 'kep', 've', 'vp','RMSE','ME']  ; parameters of extended Tofts models
	strarr_ids = ['KTRANS', 'KEP', 'VE', 'VP', 'RMSE','ME']
	sd_valid = [1,1,1,1,0,0]
	strarr_idsinf[0:2,*] = [ $
		['Ktrans =',  ' min-1',  'Ktrans (min-1)'],$
		['kep =',     ' min-1',  'kep (min-1)'   ],$
		['ve =',      ' ',       've'            ],$
		['vp =',      ' ',       'vp'            ],$
		['RMSE =',  ' mM (Ct)',  'RMSE  (mM)'  ],$
		['ME =',    ' mM (Ct)',  'ME (mM)']]
	strarr_idsinf[3,nl++] = 'Volume transfer constant between blood plasma and EES (extended Tofts model)' ; KTRANS
	strarr_idsinf[3,nl++] = 'Rate constant between EES and blood plasma (extended Tofts model)'
	strarr_idsinf[3,nl++] = 'Extracellular-extravascular instertitial space (EES) per unit volume of tissue (extended Tofts model)'
	strarr_idsinf[3,nl++] = 'Plasma volume per unit volume of tissue (extended Tofts model)'
	strarr_idsinf[3,nl++] = 'Root mean squared error (RMSE) Tofts fitting'
	strarr_idsinf[3,nl++] = 'Mean absolute error (ME) of Tofts fitting'

END
'LARSSON'   : BEGIN
	strarr_droplparams = ['S(l)', 'kep',   'RMSE','ME']             ; parameters of Larsson model
	strarr_ids = ['SL', 'KEP', 'RMSE','ME']
	sd_valid = [1,1,0,0] ; Parameters with calculation of standard deviation enabled (=1)
	strarr_idsinf[0:2,*] = [ $
		['S(l) =',      ' n.d',   'S(l) (n.d)' ],$
		['kep(l) =',    ' min-1', 'kep (min-1)'],$
		['RMSE =',  ' RCE (%)',   'RMSE (%)'], $
		['ME =',    ' RCE (%)',   'ME (%)']]

	strarr_idsinf[3,nl++] = 'Initial slope of signal (Larsson model)'
	strarr_idsinf[3,nl++] = 'Rate wash-in constant between EES and blood plasma (Larsson model)' ; kep (l)
	strarr_idsinf[3,nl++] = 'Root mean squared error (RMSE) of Larsson fitting)'
	strarr_idsinf[3,nl++] = 'Mean absolute error (ME) of Larsson fitting'

END
'RR' : BEGIN
	strarr_droplparams = ['Ktrans', 'kep', 've','RMSE','ME'] ; RR model
	strarr_ids = ['KTRANS', 'KEP', 'VE', 'RMSE','ME']
	sd_valid = [1,1,1,0,0] ; Parameters with calculation of standard deviation enabled (=1)
	strarr_idsinf[0:2,*] = [ $
		['Ktrans =',  ' min-1',  'Ktrans (min-1)'],$
		['kep =',     ' min-1',  'kep (min-1)'   ],$
		['ve =',      ' ',       've'            ],$
		['RMSE =',  ' mM (Ct)',  'RMSE  (mM)'  ],$
		['ME =',    ' mM (Ct)',  'ME (mM)']]

	strarr_idsinf[3,nl++] = 'Volume transfer constant between blood plasma and EES (RR model)' ; KTRANS
	strarr_idsinf[3,nl++] = 'Rate constant between EES and blood plasma (RR model)'
	strarr_idsinf[3,nl++] = 'Extracellular-extravascular instertitial space (EES) per unit volume of tissue (RR model)'
	strarr_idsinf[3,nl++] = 'Root mean squared error (RMSE) RR fitting'
	strarr_idsinf[3,nl++] = 'Mean absolute error (ME) of RR fitting'

END
ENDCASE

;----------------------------------------------------------------------
st_parainf  = { $
    nparam      : nparam,    $
    droplparams : strarr_droplparams, $
    ids         : strarr_ids,   $
    idsinf      : strarr_idsinf,$
    model       : model, $
    str_model   : str_model, $
    nmodel      : nmodel, $
    typemodel   : typemodel,$
    using_t10   : using_t10,$
    using_ct    : using_ct, $
    sd_valid    : sd_valid}

RETURN, st_parainf

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fStruct_Create, RESULTSTRUCT=resultstruct, ROISTRUCT=roistruct, COPY=copy, ST_IN=st_in

opt_copy=  KEYWORD_SET(copy)
alloc = opt_copy EQ 0

;-------------------------------------------------------------------------
IF KEYWORD_SET(roistruct) THEN BEGIN

	struct = { $
		opt_typeROI   : 0L, $                             ; tipo de ROI: [0, cuadrada, 1: toda la imagen, 2: libre]
		idx_points    : PTR_NEW(ALLOCATE_HEAP=alloc),$    ; indices de los puntos del ROI ( necesario en ROI free)
		idx_points_lr : PTR_NEW(ALLOCATE_HEAP=alloc),$    ; indices de los puntos en baja resolución (data_Signal_lr)
		im_idxs       : PTR_NEW(ALLOCATE_HEAP=alloc),$
		mask_lr       : PTR_NEW(ALLOCATE_HEAP=alloc),$    ; mascara 2D que incluye a la ROI
		pr_ini        : INTARR(2), $               ; Punto inicial del ROI (alta resolución)  FORMATO IDL
		roi_sz        : INTARR(2), $               ; Tamaño del rectángulo que engloba la ROI,en alta resolución
		np_lr         : 1L, $                      ; Pixeles ROI /pixel HR, la resolución de la ROI (1,2,4..)
		slicez        : 0L  $                      ; Ojo, el slice a donde va el ROI (para poder hacer ROIs en diferentes slices...)
	}
	IF opt_copy EQ 1 THEN BEGIN
		struct.opt_typeROI   = st_in.opt_typeROI
		struct.idx_points    = PTR_NEW(*st_in.idx_points)
		struct.idx_points_lr = PTR_NEW(*st_in.idx_points_lr)
		struct.im_idxs       = PTR_NEW(*st_in.im_idxs)
		struct.mask_lr       = PTR_NEW(*st_in.mask_lr)
		struct.pr_ini        = st_in.pr_ini
		struct.roi_sz        = st_in.roi_sz
		struct.np_lr         = st_in.np_lr
		struct.slicez        = st_in.slicez
	ENDIF
ENDIF

IF KEYWORD_SET(resultstruct) THEN BEGIN

	struct = { $
		dataptr        : PTRARR(7,2, ALLOCATE_HEAP=alloc),$  ; pointer to parameter  maps
		RCE_im         : PTR_NEW(ALLOCATE_HEAP=alloc),$ 	 ; pointer to RCE image  (relative constrast enhancement) (an image is a 3D data)
		arr_time       : PTR_NEW(ALLOCATE_HEAP=alloc),$      ; time array in minutes
		arr_Signal_AIF : PTR_NEW(ALLOCATE_HEAP=alloc),$
		arr_signal_RR  : PTR_NEW(ALLOCATE_HEAP=alloc),$
		arr_time_AIF   : PTR_NEW(ALLOCATE_HEAP=alloc),$
		arr_time_RR    : PTR_NEW(ALLOCATE_HEAP=alloc),$
		arr_Cp         : PTR_NEW(ALLOCATE_HEAP=alloc),$      ; Contrast agent concentration in plasma (Cp) estimated from AIF (using MRI parameters) or from analytical parameters
		data_Ct        : PTR_NEW(ALLOCATE_HEAP=alloc),$      ; Contrast agent concentration in tissue (Ct)
		data_signal_RCE: PTR_NEW(ALLOCATE_HEAP=alloc),$      ; Relative contrast enhancement used in Larsson and Hoffmann
		cp_model       : PTR_NEW(ALLOCATE_HEAP=alloc),$      ; analytical model for cp
		arr_Ct_RR      : PTR_NEW(ALLOCATE_HEAP=alloc ),$     ; Contrast agent concentration in tissue (Ct) reference region (rr)
		data_Signal_lr    : PTR_NEW(ALLOCATE_HEAP=alloc),$   ; data in ROI (low resolution, lr..)
		data_Signal_ROI   : PTR_NEW(ALLOCATE_HEAP=alloc),$   ; mean value in ROI
		data_Signal_ROI_outliers : PTR_NEW(ALLOCATE_HEAP=alloc),$   ;mean value in ROI without outliers
		scales            : FLTARR(2,7,2), $
		tofts_type     : -1}                     ; Tofft adjust -> 1: convolution. 2:  Cp biexponential

	IF opt_copy EQ 1 THEN BEGIN
		FOR i=0l, N_ELEMENTS(st_in.dataptr)-1 DO struct.dataptr[i] = PTR_NEW(*(st_in.dataptr[i]))
		struct.RCE_im     =    PTR_NEW(*st_in.RCE_im)
		struct.arr_time   =   	PTR_NEW(*st_in.arr_time)
		struct.arr_time   =   	PTR_NEW(*st_in.arr_time)
		struct.arr_Signal_AIF = PTR_NEW(*st_in.arr_Signal_AIF)
		struct.arr_signal_RR  = PTR_NEW(*st_in.arr_signal_RR)
		struct.arr_time_AIF   = PTR_NEW(*st_in.arr_time_AIF)
		struct.arr_time_RR    = PTR_NEW(*st_in.arr_time_RR)
		struct.arr_Cp      	  = PTR_NEW(*st_in.arr_Cp)
		struct.data_Ct      	=PTR_NEW(*st_in.data_Ct)
		struct.data_signal_RCE	=PTR_NEW(*st_in.data_signal_RCE)
		struct.cp_model     	=PTR_NEW(*st_in.cp_model)
		struct.arr_Ct_RR    	=PTR_NEW(*st_in.arr_Ct_RR)
		struct.data_Signal_lr   =PTR_NEW(*st_in.data_Signal_lr)
		struct.data_Signal_ROI  =PTR_NEW(*st_in.data_Signal_ROI)
		struct.data_Signal_ROI_outliers = PTR_NEW(*st_in.data_Signal_ROI_outliers)
		struct.scales         =  st_in.scales
		struct.tofts_type     =  st_in.tofts_type
	ENDIF
ENDIF
RETURN, struct

END


;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_Activate_bases, Pst, $
	ALLOW_ANALYSIS=allow_analysis,$
	FORBID_ANALYSIS=forbid_analysis,$
	OPEN_DATA=open_data,$
	OPEN_T1MAP=open_t1map,$
	FREE_DATA=free_data, $
	FREE_T1MAP=free_t1map,$
	NO_DRAW_ROI=no_draw_roi, $
	DRAW_ROI=draw_roi,$
	NEW_ROI=new_roi, $
	ROI_ANALYSIS=roi_analysis, $
	CLOSE_ROI=close_roi,$
	OPEN_ROI=open_roi,$
	IMPORT_ROI=import_roi,$
	CHANGE_ROI_TYPE=change_roi_type,$
	INITIAL_STATE=initial_state,$
	ACTIVATE_PARAMS=activate_params,$
	OPEN_MASK=open_mask,$
	FREE_MASK=free_mask


;----------------------------------------------------------------------------------------
IF KEYWORD_SET(forbid_analysis) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbase_analysis, SENSITIVE = 0
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROI, SENSITIVE= 0
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROIkinetics, SENSITIVE= 0
	(*Pst).flag_analysis = 0
	RETURN, 1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(allow_analysis) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbase_analysis, SENSITIVE=1
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROI, SENSITIVE= 1
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROIkinetics, SENSITIVE= 1
    (*Pst).flag_analysis = 1
	RETURN, 1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(change_roi_type) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbase_paramROI,   SENSITIVE= (*Pst).roi.opt_typeROI EQ 0
	WIDGET_CONTROL, (*Pst).wd.wbase_resolroi,   SENSITIVE= 1
	WIDGET_CONTROL, (*Pst).wd.wbase_analysis,   SENSITIVE= (*Pst).flag_analysis
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROI,  SENSITIVE= (*Pst).flag_analysis
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROIkinetics, SENSITIVE= (*Pst).flag_analysis
	FOR i=0, N_ELEMENTS((*Pst).wd.wbttns_viewbarROIs)-1 DO $
		WIDGET_CONTROL, (*Pst).wd.wbttns_viewbarROIs[i], SENSITIVE = (*Pst).flag_analysis
	RETURN, 1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(import_roi) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbase_paramroi,   SENSITIVE = (*Pst).roi.opt_typeROI EQ 0
	WIDGET_CONTROL, (*Pst).wd.wbase_resolroi,   SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wbase_analysis,   SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROI,   SENSITIVE= 1
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROIkinetics, SENSITIVE= 1
	FOR i=0, N_ELEMENTS((*Pst).wd.wbttns_viewbarROIs)-1 DO $
		WIDGET_CONTROL, (*Pst).wd.wbttns_viewbarROIs[i], SENSITIVE = 1
	RETURN, 1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(open_data) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbase_smain,      SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wslide_n,         SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wbase_parameters, SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wExportBar,       SENSITIVE = 1
    WIDGET_CONTROL, (*Pst).wd.woptionsbar,      SENSITIVE = 1
    WIDGET_CONTROL, (*Pst).wd.wviewbar,         SENSITIVE=1
    WIDGET_CONTROL, (*Pst).wd.wviewbar_data1,   SENSITIVE=1
    WIDGET_CONTROL, (*Pst).wd.wviewbar_data2,   SENSITIVE=1
    WIDGET_CONTROL, (*Pst).wd.wopenbar_t1free,  SENSITIVE=1
    WIDGET_CONTROL, (*Pst).wd.wopenbar_t1map,   SENSITIVE=1
	;WIDGET_CONTROL, (*Pst).wd.wopenbar_mask,   SENSITIVE=1
	WIDGET_CONTROL, (*Pst).wd.wopenbar_AIF,     SENSITIVE=1
	WIDGET_CONTROL, (*Pst).wd.wopenbar_RR,      SENSITIVE=1
  	RETURN, 1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(free_data) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbase_smain,       SENSITIVE = 0
	WIDGET_CONTROL, (*Pst).wd.wslide_n,          SENSITIVE = 0
	WIDGET_CONTROL, (*Pst).wd.wbase_parameters,  SENSITIVE = 0
	WIDGET_CONTROL, (*Pst).wd.wExportBar,        SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.woptionsbar,       SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wviewbar,          SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wviewbar_data1,    SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wviewbar_data2,    SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wopenbar_t1free,   SENSITIVE = 0
    ;WIDGET_CONTROL, (*Pst).wd.wopenbar_mask,     SENSITIVE = 0
    ;WIDGET_CONTROL, (*Pst).wd.wopenbar_maskfree, SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wopenbar_AIF,      SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wopenbar_RR,       SENSITIVE = 0

    WIDGET_CONTROL, (*Pst).wd.wbttn_drawROI, SET_BUTTON = (*Pst).flag_drawroi

   	RETURN, 1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(open_t1map) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wopenbar_t1mapfree, SENSITIVE=1
	WIDGET_CONTROL, (*Pst).wd.wviewbar_t1map,     SENSITIVE=1

	WIDGET_CONTROL, (*Pst).wdp.wtext_t10_tissue, SET_VALUE='T10 map', SENSITIVE=0
	WIDGET_CONTROL, (*Pst).wdp.wtext_t10_blood,  SET_VALUE='T10 map', SENSITIVE=0

	RETURN, 1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(free_t1map) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wopenbar_t1mapfree, SENSITIVE=0
	WIDGET_CONTROL, (*Pst).wd.wviewbar_t1map,     SENSITIVE=0

	str_T10_tissue = STRTRIM(STRING(FIX((*Pst).par.t10_tissue*1000)),2)
	str_T10_tissue = STRTRIM(STRING(FIX((*Pst).par.t10_blood*1000)),2)

	WIDGET_CONTROL, (*Pst).wdp.wtext_t10_tissue, SET_VALUE=str_T10_tissue, SENSITIVE=1
	WIDGET_CONTROL, (*Pst).wdp.wtext_t10_blood,  SET_VALUE=str_T10_blood, SENSITIVE=1

	RETURN, 1
ENDIF
 ;----------------------------------------------------------------------------------------
IF KEYWORD_SET(open_roi) THEN BEGIN
	FOR i=0l, N_ELEMENTS((*Pst).wd.wbttns_viewbarROIs)-1 DO $
		WIDGET_CONTROL, (*Pst).wd.wbttns_viewbarROIs[i], SENSITIVE = 0
	RETURN, 1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(close_roi) THEN BEGIN
	FOR i=0, N_ELEMENTS((*Pst).wd.wbttns_viewbarROIs)-1 DO $
		WIDGET_CONTROL, (*Pst).wd.wbttns_viewbarROIs[i], SENSITIVE = 1
	RETURN, 1
ENDIF
;----------------------------------------------------------------------------------------
;IF KEYWORD_SET(open_mask) THEN BEGIN
;	WIDGET_CONTROL, (*Pst).wd.wopenbar_maskfree, SENSITIVE=1
;	RETURN, 1
;ENDIF
;----------------------------------------------------------------------------------------
;IF KEYWORD_SET(free_mask) THEN BEGIN
;	WIDGET_CONTROL, (*Pst).wd.wopenbar_maskfree, SENSITIVE=0
;	RETURN, 1
;ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(no_draw_roi) THEN BEGIN
	(*Pst).vi1.Opolyline_ROI -> SetProperty,    HIDE=1
    (*Pst).vi2.Opolyline_ROI -> SetProperty,    HIDE=1
    (*Pst).vi1.Opolyline_AUX -> SetProperty,    HIDE=1
    (*Pst).vi2.Opolyline_AUX -> SetProperty,    HIDE=1
    WIDGET_CONTROL, (*Pst).wd.wbttn_newROI,    SENSITIVE=0
    WIDGET_CONTROL, (*Pst).wd.wbase_typeROI,   SENSITIVE=0
	WIDGET_CONTROL, (*Pst).wd.wbase_paramROI,  SENSITIVE=0
	WIDGET_CONTROL, (*Pst).wd.wbase_resolROI,  SENSITIVE=0
	WIDGET_CONTROL, (*Pst).wd.wbase_analysis,  SENSITIVE=0
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROI, SENSITIVE=0
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROIkinetics, SENSITIVE= 0
	WIDGET_CONTROL, (*Pst).wd.wtext_info2,    SET_VALUE=''

	WIDGET_CONTROL, (*Pst).wd.wslide_z, SENSITIVE = 1
	FOR i=0, N_ELEMENTS((*Pst).wd.wbttns_viewbarROIs)-1 DO $
		WIDGET_CONTROL, (*Pst).wd.wbttns_viewbarROIs[i], SENSITIVE = 0

	RETURN,1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(draw_roi) THEN BEGIN
	(*Pst).vi1.Opolyline_ROI -> SetProperty, HIDE=0
	(*Pst).vi2.Opolyline_ROI -> SetProperty, HIDE=0
	(*Pst).vi1.Opolyline_AUX -> SetProperty, HIDE=0
    (*Pst).vi2.Opolyline_AUX -> SetProperty, HIDE=0
    WIDGET_CONTROL, (*Pst).wd.wbttn_newROI,     SENSITIVE=1
	WIDGET_CONTROL, (*Pst).wd.wbase_typeROI,    SENSITIVE=1
	WIDGET_CONTROL, (*Pst).wd.wbase_paramROI,   SENSITIVE=(*Pst).roi.opt_typeROI EQ 0
	WIDGET_CONTROL, (*Pst).wd.wbase_resolROI,   SENSITIVE=1
	WIDGET_CONTROL, (*Pst).wd.wbase_analysis,   SENSITIVE=(*Pst).flag_analysis
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROI,  SENSITIVE=(*Pst).flag_analysis
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROIkinetics, SENSITIVE= (*Pst).flag_analysis

	FOR i=0, N_ELEMENTS((*Pst).wd.wbttns_viewbarROIs)-1 DO $
		WIDGET_CONTROL, (*Pst).wd.wbttns_viewbarROIs[i], SENSITIVE = (*Pst).flag_freeROIclose

	RETURN,1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(roi_analysis) THEN BEGIN

	FOR i=0l, N_ELEMENTS((*Pst).wd.wbases_params[*,1])-1 DO BEGIN
		WIDGET_CONTROL, (*Pst).wd.wbases_params[i,0], SENSITIVE = 0
		WIDGET_CONTROL, (*Pst).wd.wbases_params[i,1], SENSITIVE = (*Pst).flags_applyparams[i]
		WIDGET_CONTROL, (*Pst).wd.wbases_params[i,2], SENSITIVE = 0
	ENDFOR

	WIDGET_CONTROL, (*Pst).wd.wbase_typeroi,       SENSITIVE=  0
	WIDGET_CONTROL, (*Pst).wd.wbase_paramroi,      SENSITIVE=  0
	WIDGET_CONTROL, (*Pst).wd.wbase_resolroi,      SENSITIVE=  1
	WIDGET_CONTROL, (*Pst).wd.wbttn_resetROIresult, SENSITIVE= 1
	WIDGET_CONTROL, (*Pst).wd.wbarmenu_limitframes,SENSITIVE=  0
	WIDGET_CONTROL, (*Pst).wd.wexportbar,          SENSITIVE=  1
	WIDGET_CONTROL, (*Pst).wd.wbttn_hideroilimits, SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wslide_z, SENSITIVE = 0
	FOR i=0, N_ELEMENTS((*Pst).wd.wbttns_viewbarROIs)-1 DO $
		WIDGET_CONTROL, (*Pst).wd.wbttns_viewbarROIs[i], SENSITIVE = 1

	WIDGET_CONTROL, (*Pst).wd.wbttn_drawROI, SENSITIVE = 0; NEW february 2012

	RETURN,1

ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(no_image) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbase_smain,   SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wExportBar,    SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.woptionsbar,   SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wviewbar,      SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wviewbar_data1, SENSITIVE = 0
    WIDGET_CONTROL, (*Pst).wd.wviewbar_data2, SENSITIVE = 0
       FOR i=0, N_ELEMENTS((*Pst).wd.wbttns_viewbarROIs)-1 DO $
		WIDGET_CONTROL, (*Pst).wd.wbttns_viewbarROIs[i], SENSITIVE = 0
    RETURN,1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(new_roi) THEN BEGIN

	WIDGET_CONTROL, (*Pst).wd.wbase_typeroi,    SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wbase_paramroi,   SENSITIVE = (*Pst).roi.opt_typeROI EQ 0
	WIDGET_CONTROL, (*Pst).wd.wbase_resolroi,   SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wbase_analysis,   SENSITIVE = 0
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROI,  SENSITIVE=  0
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROIkinetics, SENSITIVE= 0
	WIDGET_CONTROL, (*Pst).wd.wslide_z, SENSITIVE    = 1

	FOR i=0, N_ELEMENTS((*Pst).wd.wbttns_viewbarROIs)-1 DO $
		WIDGET_CONTROL, (*Pst).wd.wbttns_viewbarROIs[i], SENSITIVE = 0

	RETURN,1
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(initial_state) THEN BEGIN

	; state when interface begins

	WIDGET_CONTROL, (*Pst).wd.wbase_smain,      SENSITIVE = 0 ; Base que agrupa las acciones sobre una imagen
	WIDGET_CONTROL, (*Pst).wd.wbase_analysis,   SENSITIVE = 0 ; Base que agrupa el analisis de la roi
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROI,  SENSITIVE = 0
	WIDGET_CONTROL, (*Pst).wd.wbttn_exportROIkinetics, SENSITIVE= 0
	WIDGET_CONTROL, (*Pst).wd.wbase_parameters, SENSITIVE = 0 ; Base con los tabs de parámetros
	WIDGET_CONTROL, (*Pst).wd.wopenbar_t1map,     SENSITIVE = 0
	WIDGET_CONTROL, (*Pst).wd.wdropl_typeROI,   SET_DROPLIST_SELECT=(*Pst).roi.opt_typeROI
	WIDGET_CONTROL, (*Pst).wd.wbase_paramROI,   SENSITIVE= (*Pst).roi.opt_typeROI EQ 0
	WIDGET_CONTROL, (*Pst).wd.wbase_resolroi,   SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wbase_newROI,     SENSITIVE = 1
	WIDGET_CONTROL, (*Pst).wd.wbttn_interp,     SET_BUTTON = (*Pst).flag_interp
	WIDGET_CONTROL, (*Pst).wd.wopenbar_AIF,     SENSITIVE = 0
	WIDGET_CONTROL, (*Pst).wd.wopenbar_RR,      SENSITIVE = 0

	pos_resolution = (WHERE((*Pst).wd.wdropl_ROIresolution_Array EQ (*Pst).roi.np_lr))[0]
	WIDGET_CONTROL, (*Pst).wd.wdropl_ROIresolution, SET_DROPLIST_SELECT=pos_resolution
	;----------------------------------------------------------------------
	WIDGET_CONTROL, (*Pst).wd.wdropls_model, SET_DROPLIST_SELECT=0
	;----------------------------------------------------------------------
	WIDGET_CONTROL, (*Pst).wd.wbttn_drawROI, SET_BUTTON = (*Pst).flag_drawroi
ENDIF
;----------------------------------------------------------------------------------------
IF KEYWORD_SET(activate_params) THEN BEGIN

	FOR i=0l, N_ELEMENTS((*Pst).wd.wbases_params[*,1])-1 DO BEGIN
		WIDGET_CONTROL, (*Pst).wd.wbases_params[i,1], SENSITIVE=(*Pst).flags_applyparams[i]
		WIDGET_CONTROL, (*Pst).wd.wbases_params[i,2], SENSITIVE=(*Pst).flags_applyparams[i]
	ENDFOR
ENDIF

RETURN,1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_SuffixResultName, ST_ROI=st_roi, SIZEIMAGE=sizeimage, str_add

IF st_roi.opt_typeROI NE 2 THEN BEGIN
	point_v2 = get_point_YupConvention(st_roi.pr_ini, DIMENSIONS=sizeimage, SIZE_ROI=st_roi.roi_sz)
	str_pointini = '.pi[' + STRTRIM(point_v2[0],2) + ',' + STRTRIM(point_v2[1],2) + ']'
	str_pointini+= '.s[' + STRTRIM(st_roi.roi_sz[0]/st_roi.np_lr[0],2) + ',' + $
		STRTRIM(st_roi.roi_sz[1]/st_roi.np_lr[0],2) + ']'
	;str_pointini+= '.r[' + STRTRIM((*Pst).roi.np_lr[0],2)  + ',' + STRTRIM((*Pst).roi.np_lr[0],2)  + ']'
	str_pointini+= '.r[' + STRTRIM(st_roi.np_lr[0],2)  + ']'
ENDIF ELSE BEGIN
	str_pointini = str_add
	;str_pointini+= '.r[' + STRTRIM((*Pst).results.np_lr[0],2)  + ',' + STRTRIM((*Pst).results.np_lr[0],2)  + ']'
	str_pointini+= '.r[' + STRTRIM(st_roi.np_lr[0],2)  + ']'
ENDELSE

RETURN, str_pointini

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_BoxROI_set, Pst, POINT_INI=point_ini, DRAWTWO=drawtwo, ID_TEXT=id_text

	;point_ini IDL standard ("y-down"): esquina inferior izquierda
	;"Y-up" criteria is only to represent in screen ...

	n_ppx = (*Pst).inf.size_roi[0]
	n_ppy = (*Pst).inf.size_roi[1]
	tam = [(*Pst).inf.size_im[0],(*Pst).inf.size_im[1],(*Pst).inf.nslices,(*Pst).inf.nframes]


	IF ((point_ini[0]) LT 0) OR ((point_ini[1]) LT 0) OR $
	   ((point_ini[0]) GE tam[0]) OR ((point_ini[1]) GE tam[1]) THEN RETURN,-1

	(*Pst).inf.p_ini = point_ini
	;---------------------------------------------------------------
	pl_1x = point_ini[0]
	pl_1y = point_ini[1]

	pl_2x = point_ini[0]+n_ppx
	pl_2y = point_ini[1]

	pl_3x = point_ini[0]+n_ppx
	pl_3y = point_ini[1]+n_ppy

	pl_4x = point_ini[0]
	pl_4y = point_ini[1]+n_ppy

	data_polyline =[[pl_1x,pl_1y],[pl_2x,pl_2y],[pl_3x,pl_3y],[pl_4x,pl_4y],[pl_1x,pl_1y]]
	PRINT,''
	PRINT,'[' + STRTRIM(pl_1x,2) + ', ' + STRTRIM(pl_1y,2) + ']' + '  ' + $
		'[' + STRTRIM(pl_2x,2) + ', ' + STRTRIM(pl_2y,2) + ']'  +  '  ' + $
		'[' + STRTRIM(pl_3x,2) + ', ' + STRTRIM(pl_3y,2) + ']'  +  '  ' + $
		'[' + STRTRIM(pl_4x,2) + ', ' + STRTRIM(pl_4y,2) + ']'

	(*Pst).vi2.Opolyline_ROI -> SetProperty, DATA=data_polyline
	;--------------------------------------------------------------


	p_imf_1 = get_point_YupConvention([pl_4x,pl_4y-1], DIMENSIONS=(*Pst).inf.size_im)
	p_imf_2 = get_point_YupConvention([pl_2x-1,pl_2y], DIMENSIONS=(*Pst).inf.size_im)

	str_text = '[' +  STRTRIM(p_imf_1[0],2) + ', ' + STRTRIM(p_imf_1[1],2) + $
		'] -> [' + STRTRIM(p_imf_2[0],2) + ', ' + STRTRIM(p_imf_2[1],2) + ']'

	FOR i=0l, N_ELEMENTS(id_text)-1 DO BEGIN
		WIDGET_CONTROL, id_text[i], SET_VALUE=str_text
	ENDFOR

	WIDGET_CONTROL, (*Pst).wd.wtext_roi_px, SET_VALUE=STRTRIM(STRING(p_imf_1[0]),2)
	WIDGET_CONTROL, (*Pst).wd.wtext_roi_py, SET_VALUE=STRTRIM(STRING(p_imf_1[1]),2)

	IF KEYWORD_SET(drawtwo) THEN BEGIN
		(*Pst).vi1.Opolyline_ROI -> SetProperty, DATA=data_polyline

	ENDIF

	RETURN, (*Pst).inf.p_ini

END

;**************************************************************************************************
;**************************************************************************************************

PRO fwdcemri_Change_palette, DATA=data

    TVLCT, R,G,B, /GET

    CASE data.changing_palette OF
    1 : BEGIN
	    data.vi1.opalette   -> setproperty, RED=r,GREEN=G, BLUE=b
	    data.vi1.ohcolorbar -> setproperty, PALETTE  = data.vi1.opalette
	    WIDGET_CONTROL, data.wd.wDraw1, GET_VALUE = window_1
    	window_1 -> Draw, data.vi1.oView
	END
	2 : BEGIN
	   	data.vi2.opalette   -> setproperty, RED=r,GREEN=G, BLUE=b
    	data.vi2.ohcolorbar -> setproperty, PALETTE  = data.vi2.opalette
      	WIDGET_CONTROL, data.wd.wDraw2, GET_VALUE = window_2
    	window_2 -> Draw, data.vi2.oView
	END
	ELSE:
    ENDCASE

 END

;****************************************************************************
;****************************************************************************

FUNCTION fwdcemri_CheckDrawRoi, Pst

	IF (N_ELEMENTS(*(*Pst).data) GT 2) THEN BEGIN

	ENDIF
	RETURN, 1
END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_DataFromSubidx, ST_ROI=st_roi, image

min_p = st_roi.pr_ini
max_p = min_p + st_roi.roi_sz-1
subim = image[min_p[0]:max_p[0],min_p[1]:max_p[1]]

IF st_roi.np_lr GT 1 THEN BEGIN
	dim_subim = (max_p-min_p+1)/st_roi.np_lr
	subim_lr = REBIN(subim, dim_subim[0], dim_subim[1])
ENDIF ELSE BEGIN
	subim_lr = subim
ENDELSE

data_Signal_lr = subim_lr[*st_roi.idx_points_lr]
RETURN, data_Signal_lr

END

;****************************************************************************
;****************************************************************************

FUNCTION fwdcemri_Draw_Image, Pst, CONTROL_ZOOM=control_zoom, ERROR_STR = error_str, NODRAW=nodraw

	opt_nodraw = KEYWORD_SET(nodraw) ; this option not redraw the window. Useful for "simultaneous redraw of many windows)

     error_str = 'Error in module "fwdcemri_Draw_image"'

    IF N_PARAMS() NE 1 THEN RETURN, -1
    opt_controlzoom = KEYWORD_SET(control_zoom)

    ;------------------------------------------------
    WIDGET_CONTROL,(*Pst).wd.wslide_z,  GET_VALUE=n_z
    WIDGET_CONTROL,(*Pst).wd.wslide_n,  GET_VALUE=n_n
    (*Pst).inf.current_slice = n_z
    (*Pst).inf.current_frame = n_n
    WIDGET_CONTROL,(*Pst).wd.wlabel_z, SET_VALUE=STRTRIM(n_z+1,2)
    WIDGET_CONTROL,(*Pst).wd.wlabel_n, SET_VALUE=STRTRIM(n_n+1,2)

	tam = SIZE(*(*Pst).data, /DIMENSIONS)
    img = REFORM((*(*Pst).data)[*,*,n_z,n_n])


    (*Pst).vi1.omodel  -> setProperty, HIDE=0
    (*Pst).vi1.oimage  -> setProperty, INTERPOLATE=(*Pst).flag_interp

	min_val =  MIN(img,    MAX=max_val)
    img_scl =  BYTSCL(img, MIN=min_val, MAX=max_val)

    xMax = tam[0] & yMax = tam[1]

    (*Pst).vi1.oimage -> setProperty, DATA = img_scl, DIMENSIONS=[xMax, yMax], GREYSCALE=0,$
    	SUB_RECT = [0,0,xMax, yMax], HIDE=0

    ;---------------------------------------------------
    IF opt_controlzoom EQ 1 THEN BEGIN
    	rel_xy = tam[0]*1d/(tam[1] + tam[0]*(*Pst).opt.colorbar_relativesize[0])
    	fwidget_ControlZoomDraw, (*Pst).wd.wdraw1, REL_XY=rel_xy, ZOOM=control_zoom, XLIMITS=(*Pst).opt.view_size_limits_x, YLIMITS=(*Pst).opt.view_size_limits_Y, SCROLL_LIMIT=(*Pst).opt.view_scroll_limit
    ENDIF
    ;---------------------------------------------------

    IF (opt_nodraw EQ 0) THEN BEGIN
    	WIDGET_CONTROL, (*Pst).wd.wDraw1, GET_VALUE = window_1
    	window_1 -> Draw, (*Pst).vi1.oView
    ENDIF

    undefine, error_str
	RETURN , 1

END

;****************************************************************************
;****************************************************************************

FUNCTION fwdcemri_Draw_RCEImage, Pst, CONTROL_ZOOM=control_zoom, ERROR_STR = error_str, $
	CALCULATE=calculate, NODRAW=nodraw

	opt_nodraw = KEYWORD_SET(nodraw) ; this option not redraw the window. Useful for "simultaneous redraw of many windows)

    error_str = 'Error in module "fwdcemri_Draw_RCEimage"'

    IF N_PARAMS() NE 1 THEN RETURN, -1
    IF KEYWORD_SET(control_zoom) THEN opt_controlzoom=1 ELSE opt_controlzoom=0

	WIDGET_CONTROL, /HOURGLASS
    ;--------------------------------------------------
    WIDGET_CONTROL,(*Pst).wd.wslide_z,  GET_VALUE=n_z
    WIDGET_CONTROL,(*Pst).wd.wslide_n,  GET_VALUE=n_n
    (*Pst).inf.current_slice = n_z
    (*Pst).inf.current_frame = n_n
    WIDGET_CONTROL,(*Pst).wd.wlabel_z, SET_VALUE=STRTRIM(STRING(n_z+1),2)
    WIDGET_CONTROL,(*Pst).wd.wlabel_n, SET_VALUE=STRTRIM(STRING(n_n+1),2)

	WIDGET_CONTROL,(*Pst).wd.wslide_minval, GET_VALUE=minval_percent
    WIDGET_CONTROL,(*Pst).wd.wslide_maxval, GET_VALUE=maxval_percent
    maxval_percent = maxval_percent > (minval_percent+1)
    WIDGET_CONTROL,(*Pst).wd.wslide_maxval, SET_VALUE=maxval_percent

    WIDGET_CONTROL,(*Pst).wd.wslide_threshold, GET_VALUE=threshold_percent

    WIDGET_CONTROL,(*Pst).wd.wtext_minval,    SET_VALUE=STRTRIM(STRING(minval_percent),2) + '%'
    WIDGET_CONTROL,(*Pst).wd.wtext_maxval,    SET_VALUE=STRTRIM(STRING(maxval_percent),2) + '%'
    WIDGET_CONTROL,(*Pst).wd.wlabel_threshold, SET_VALUE=STRTRIM(STRING(threshold_percent),2) + '%'

	tam = SIZE(*(*Pst).data, /DIMENSIONS)
    ;---------------------------------------------------

    IF KEYWORD_SET(calculate) THEN BEGIN
    	*(*Pst).results.RCE_im = FLTARR((*Pst).inf.size_im[0],(*Pst).inf.size_im[1],(*Pst).inf.nslices)
    	FOR k=0l, (*Pst).inf.nslices-1 DO BEGIN
	       	min_threshold  = MAX(REFORM((*(*Pst).data)[*,*,k,*]))*threshold_percent/100.0
	    	;max_threshold = MAX(REFORM((*(*Pst).data)[*,*,k,*]))*maxval*1.0/100
	       	;---------------------------------------------------------------------
		    img_RCE = RelativeContrastEnhancement(REFORM((*(*Pst).data)[*,*,k,*]),FRAME_INJECTION=(*Pst).par.nframe_injection, $
		    	MIN_THRESHOLD=min_threshold, MAX_THRESHOLD=max_threshold, ABSOLUTEMINMAX=(*Pst).opt_typeRCE EQ 2)
		    ;---------------------------------------------------------------------
		   	; se guarda un valor que difiere del mostrado en pantalla porque no está limitado
	    	; por minval_percent y maxval_percent
	    	(*(*Pst).results.RCE_im)[*,*,k] = img_RCE
	    ENDFOR
    ENDIF

  	;imgview_RCE =   (minval_percent/100.0) > (*(*Pst).results.RCE_im)[*,*,n_z]  < (maxval_percent/100.0)
  	;img_scl      = BYTSCL(imgview_RCE, MIN=minval_percent/100.0, MAX=maxval_percent/100.0)
  	imgview_RCE =   (*(*Pst).results.RCE_im)[*,*,n_z]
   	img_scl     =  BYTSCL(imgview_RCE)

  	;---------------------------------------------------
    IF opt_controlzoom EQ 1 THEN BEGIN
       rel_xy = tam[0]*1d/(tam[1] + tam[0]*(*Pst).opt.colorbar_relativesize[0])
       fwidget_ControlZoomDraw, (*Pst).wd.wdraw2, REL_XY=rel_xy, ZOOM=control_zoom, XLIMITS=(*Pst).opt.view_size_limits_x, YLIMITS=(*Pst).opt.view_size_limits_Y,SCROLL_LIMIT=(*Pst).opt.view_scroll_limit
    ENDIF
    ;---------------------------------------------------
    st_drawbase = WIDGET_INFO((*Pst).wd.wdraw2, /GEOMETRY)

    (*Pst).vi2.omodel  -> setProperty, HIDE=0
    (*Pst).vi2.oimage  -> setProperty, INTERPOLATE=(*Pst).flag_interp

	img_scl = BYTSCL(imgview_RCE, MIN=minval_percent/100.0, MAX=maxval_percent/100.0)

    (*Pst).vi2.ohcolorbar->setproperty, RANGE=[minval_percent,maxval_percent]

   	xMax = tam[0] & yMax = tam[1]
   	(*Pst).vi2.oimage -> setProperty, DATA = img_scl, SUB_RECT = [0,0,xMax, yMax], DIMENSIONS=[xMax, yMax], GREYSCALE=0, HIDE=0


	IF (opt_nodraw EQ 0) THEN BEGIN
		WIDGET_CONTROL, (*Pst).wd.wDraw2, GET_VALUE = window
		window -> Draw, (*Pst).vi2.oView
	ENDIF
    undefine, error_str

    RETURN , 1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_DynDataFromIdx, Pst, IDX_POINTS=idx_points

; Extract array of dynamic data from pixel indexes

IF N_ELEMENTS(idx_points) EQ 0 THEN BEGIN
	IF N_ELEMENTS(*(*Pst).roi.idx_points) LE 0 THEN BEGIN
		RETURN,-1
	ENDIF
	idxp = *(*Pst).roi.idx_points
ENDIF ELSE BEGIN
	idxp = idx_points
ENDELSE
n_points = N_ELEMENTS(idxp)

WIDGET_CONTROL,(*Pst).wd.wslide_z,  GET_VALUE=n_z
(*Pst).inf.current_slice = n_z

arr_result = FLTARR(n_points, (*Pst).inf.nframes)

FOR i=0l, (*Pst).inf.nframes-1 DO BEGIN
	data_temp = REFORM((*(*Pst).data)[*,*,n_z,i])
	arr_result[*,i] = data_temp[idxp]
ENDFOR

RETURN, arr_result

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_ExportImage, image_view, NAME=name, OPTION_NAME=option_name, PATH_RESULT=path_result, INFO=info, NO_VIEW=no_view

	IF SIZE(image_view, /N_DIMENSIONS) NE 3 THEN BEGIN
		ok = DIALOG_MESSAGE('Error saving image', /CENTER, /INFO)
		RETURN,-1
	ENDIF

    opt_win  = KEYWORD_SET(no_view) EQ 0
    opt_name = KEYWORD_SET(option_name) ; 1 no salva en disco, solo para poner la imagen y exportar a PPT

	IF opt_win EQ 1 THEN BEGIN

		xsize = (SIZE(image_view, /DIMENSIONS))[1]
		ysize = (SIZE(image_view, /DIMENSIONS))[2]

		IF opt_name EQ 1 THEN BEGIN
			ysize_new = CEIL(ysize*1.10)
			imag_new  = BYTARR(3,xsize, ysize_new)+255
			imag_new[*,*,ysize_new-ysize:ysize_new-1] = image_view
			image_view = TEMPORARY(imag_new)
			ysize = TEMPORARY(ysize_new)
		ENDIF

    	win = (!D.WINDOW+1) MOD 32
    	WINDOW, win, XSIZE=xsize, YSIZE=ysize, XPOS=20, YPOS=20
    	TV, image_view, TRUE=1

		IF opt_name EQ 1 THEN BEGIN

			DEVICE, SET_FONT = "Arial*15"

			Width  = 1.5*(!D.y_ch_size*1.0/!D.y_size) ; Espaciado y medio
			lenght = 1.0*(!D.x_ch_size)/!D.x_size*MAX(STRLEN(info))

			XYOUTS, 0.03, 0.01,       info[1], /NORMAL, FONT=0, COLOR=0
			XYOUTS, 0.03, 0.01+Width, info[0], /NORMAL, FONT=0, COLOR=0
			RETURN,-1
		ENDIF
    ENDIF

    IF N_ELEMENTS(name) NE 0 THEN BEGIN
		str_file_ini = name
	ENDIF ELSE str_file_ini = ''

	str_file = DIALOG_PICkFILE(PATH=path_result, TITLE='save TIFF file', $
		FILE=str_file_ini, /WRITE)
	IF str_file EQ '' THEN RETURN,-1
	image_view = REVERSE(image_view,3)
	WRITE_TIFF, str_file + '.tif', image_view, ORIENTATION=1
	path_result = get_path(str_file)

RETURN, image_view

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_ExportROI, st_roi, IMAGE_SIZE=image_size, PATH_RESULT=path_result, OPOLYLINE=opolyline

	reverse_lines = 1l ; be careful, must be 1 in all cases
	idx_points    = *st_roi.idx_points
    str_file_ini  = 'ROI_positions'
    opt_write_ascii = 1
    opt_write_image = 1 AND opt_write_ascii
    opt_name  = 1

	str_pointini = fwdcemri_suffixresultname(ST_ROI=st_roi, SIZEIMAGE=image_size, '[ROIf]')

	str_file_ini+=str_pointini + '.sav'

	str_file = DIALOG_PICKFILE(PATH=path_result, TITLE='Export SAV file', $
		DEFAULT_EXTENSION = 'sav', WRITE=1, FILE=str_file_ini)
	IF str_file EQ '' THEN RETURN,-1

	path_result  = get_path(str_file)
	str_file_sav = get_path(str_file) + get_name_field(str_file) + '.sav'

	IF opt_write_ascii THEN BEGIN

		idx_points_corners = IdxPoints_FirstFromlowRes(MASK_IDX=*st_roi.im_idxs, DIMENSIONS=image_size, $
			ROI_RESOLUTION=st_roi.np_lr, MASK_OUTPUT=mask_output)
		IF idx_points_corners[0] EQ -1 THEN BEGIN
			ok = DIALOG_MESSAGE('Error saving ROI, please do not change ROI resolution or size since last analysis', /INFO)
			RETURN,-1
		ENDIF

		idx_points_Yup = IdxPoints_ChangeYconvention(idx_points_corners, DIMENSIONS=image_size, ARRAY_XY=array_xy)
		str_array = Translate_data2strings(array_xy, N_DECIMALS=5, REVERSE_LINES=reverse_lines, INTEGER=1)

		; difference: in SAV, idx_points is in IDL convention, in the ascii code, is in IMAFEN convention
		file_name_txt  = get_path(str_file) + get_name_field(str_file) + '.txt'
		str_info = STRARR(4)
		str_info[0] = '## Slice size (x,y): ' + STRTRIM(image_size[0],2) +  ' ' + STRTRIM(image_size[1],2)
		str_info[1] = '## ROI resolution  : ' + STRTRIM(st_roi.np_lr[0],2)
		str_info[2] = '## x - y'
		str_array = [str_info,str_array]
		ok = WriteFile_ASCII(file_name_txt, ARRAY=str_array, ERROR_ID=error_id)
	ENDIF

	IF opt_write_image THEN BEGIN
		mask = (mask_output*127b) > 32b
		file_name_tiff = get_path(str_file) + get_name_field(str_file) + '.png'
		WRITE_PNG, file_name_tiff, mask
	ENDIF

	Opolyline->getproperty, DATA=polylineROI

	ResolutionRoi = st_roi.np_lr
	;typeROI      = (*Pst).roi.opt_typeROI
	typeROI       = 2; se escribirá como roi free, aunque venga de uno cuadrado
	;----------------------------------------------
	SAVE, idx_points, polylineROI, resolutionROI, typeROI, FILENAME=str_file_sav
	;----------------------------------------------

	IF FILE_TEST(str_file_sav) THEN BEGIN
		ok = DIALOG_MESSAGE(['Result saved in file ', get_name(str_file_sav)], /INFO)
	ENDIF
	RETURN, 1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_FullROI_set, Pst, DRAWTWO=drawtwo, ID_TEXT=id_text

		;point_ini SIGUE EL CRITERIO IDL: esquina inferior derecha
	;el criterio "Y-up" solo para respresentar en pantalla...

	n_ppx = (*Pst).inf.size_im[0]
	n_ppy = (*Pst).inf.size_im[1]
	;---------------------------------------------------------------
	; Ahora exclusivo para la polyline...
	pl_1x = 0l
	pl_1y = 0l
	pl_2x = n_ppx
	pl_2y = 0l
	pl_3x = n_ppx
	pl_3y = n_ppy
	pl_4x = 0l
	pl_4y = n_ppy
	data_polyline =[[pl_1x,pl_1y],[pl_2x,pl_2y],[pl_3x,pl_3y],[pl_4x,pl_4y],[pl_1x,pl_1y]]

	(*Pst).vi2.Opolyline_ROI -> SetProperty, DATA=data_polyline
	;--------------------------------------------------------------
	p_imf_1 = get_point_YupConvention([pl_4x,pl_4y-1], DIMENSIONS=(*Pst).inf.size_im)
	p_imf_2 = get_point_YupConvention([pl_2x-1,pl_2y], DIMENSIONS=(*Pst).inf.size_im)

	str_text = '[' +  STRTRIM(p_imf_1[0],2) + ', ' + STRTRIM(p_imf_1[1],2) + $
		'] -> [' + STRTRIM(p_imf_2[0],2) + ', ' + STRTRIM(p_imf_2[1],2) + ']'

	FOR i=0l, N_ELEMENTS(id_text)-1 DO BEGIN
		WIDGET_CONTROL, id_text[i], SET_VALUE=str_text
	ENDFOR

	IF KEYWORD_SET(drawtwo) THEN BEGIN
		(*Pst).vi1.Opolyline_ROI -> SetProperty, DATA=data_polyline
	ENDIF

	RETURN, [0,0l]


END


;****************************************************************************
;****************************************************************************

FUNCTION fwdcemri_Hide_Images, Pst

	(*Pst).vi1.oView->SetProperty,  HIDE = 1
	(*Pst).vi2.oView->SetProperty,  HIDE = 1
	(*Pst).vi1.omodel->SetProperty, HIDE = 1
	(*Pst).vi2.omodel->SetProperty, HIDE = 1
	(*Pst).vi1.oimage->SetProperty, HIDE = 1
	(*Pst).vi2.oimage->SetProperty, HIDE = 1

	(*Pst).vi1.opolyline_ROI->setproperty, HIDE = 1
	(*Pst).vi2.opolyline_ROI->setproperty, HIDE = 1

	WIDGET_CONTROL, (*Pst).wd.wDraw1, GET_VALUE = window_draw
   	window_draw -> Draw, (*Pst).vi1.oView

	WIDGET_CONTROL, (*Pst).wd.wDraw2,  GET_VALUE = window_draw
	window_draw -> Draw, (*Pst).vi2.oView

END


;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_Load_image, Pst, data, FILENAME=filename, $
         ERROR_STR = error_str, MASK=mask

    error_str = 'Error in module "fwdcemri_Load_image"'

    IF N_PARAMS() NE 2 THEN RETURN, -1

    type_format = SIZE(data, /TYPE)
    IF  SIZE(Data, /N_DIMENSIONS) NE 4 THEN RETURN, -2

    IF N_ELEMENTS(filename) THEN BEGIN
    	str_path = FILE_DIRNAME(filename, /MARK_DIRECTORY)
    	str_name = FILE_BASENAME(filename)
    ENDIF ELSE BEGIN
    	str_path = '.\'
    	str_name = 'temp'
	ENDELSE

    (*Pst).inf.name    = str_name
    (*Pst).inf.path    = str_path

	(*Pst).inf.size_im = (SIZE(data,/DIMENSIONS))[0:1]
	(*Pst).inf.nslices = (SIZE(data,/DIMENSIONS))[2]
	(*Pst).inf.nframes = (SIZE(data,/DIMENSIONS))[3]

    (*Pst).par.nframes = (*Pst).inf.nframes
    (*Pst).par.nframes_adjust = (*Pst).inf.nframes-(*Pst).par.nframe_injection

    WIDGET_CONTROL, (*Pst).wdp.wtext_nframes, SET_VALUE= STRTRIM(STRING(FIX((*Pst).par.nframes)),2)

	WIDGET_CONTROL,(*Pst).wd.wslide_z, SET_VALUE=(*Pst).inf.nslices/2, SET_SLIDER_MAX=(*Pst).inf.nslices-1, SENSITIVE=1
    WIDGET_CONTROL,(*Pst).wd.wslide_n, SET_VALUE=0,        SET_SLIDER_MAX=(*Pst).inf.nframes-1, SENSITIVE=1
    (*Pst).inf.current_slice = (*Pst).inf.nslices/2
    (*Pst).inf.current_frame = 0


	*(*Pst).results.arr_time = INDGEN((*Pst).par.nframes)*(*Pst).par.frame_period
	  ; tiempo total en minutos

    *(*Pst).data = data

	IF N_ELEMENTS(mask) NE 0 THEN BEGIN
		size_mask = SIZE(mask, /DIMENSIONS)
		IF (*Pst).inf.size_im[0]   NE size_mask[0] THEN RETURN, -1
		IF (*Pst).inf.size_im[1]   NE size_mask[1] THEN RETURN, -1
		IF (*Pst).inf.nslices NE size_mask[2] THEN RETURN, -1
		*(*Pst).mask = mask
	ENDIF ELSE BEGIN
		IF  N_ELEMENTS(*(*Pst).mask) EQ 0 THEN BEGIN
			*(*Pst).mask = BYTARR((*Pst).inf.size_im[0],(*Pst).inf.size_im[1],(*Pst).inf.nslices)
		ENDIF
	ENDELSE

	max_value = MAX(data)
	(*Pst).vi1.ohcolorbar-> setproperty, PALETTE=(*Pst).vi1.opalette, RANGE=[0,max_value], TITLE=''

    RETURN, 1
END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_ImageInfo, Pst, TITLE=title

	strarr_info=STRARR(1,8)
	n=1l
	strarr_info[n++] = 'File name:  ' + (*Pst).inf.name
	strarr_info[n++] = 'Path:       ' + (*Pst).inf.path
	strarr_info[n++] = 'Size (x-y): ' + Translate_numberarray2string(ARRAY=(*Pst).inf.size_im, SEP=',', BRACKETS=1) + ' pixels'
	strarr_info[n++] = 'Slices (z): ' + STRTRIM(STRING((*Pst).inf.nslices),2)
	strarr_info[n++] = 'temporal frames: ' + STRTRIM(STRING((*Pst).inf.nframes),2)
	n++
	strarr_info[n++] = 'Frame period: ' + STRTRIM(STRING((*Pst).par.frame_period, FORMAT='(F10.4)'),2) + ' min. (' + $
									STRTRIM(STRING((*Pst).par.frame_period*60, FORMAT='(F10.4)'),2) + ' s)'

	max_length= (MAX(STRLEN(strarr_info))+50) < 200
	title = 'Information' + STRJOIN(REPLICATE('_', max_length))

	RETURN, strarr_info


END

;**************************************************************************************************
;**************************************************************************************************


FUNCTION fwdcemri_PlotsOpenedType, id, wbttn_plotsOpened, opt_plotsOpened


	IF (opt_plotsOpened EQ 1) AND (id EQ wbttn_plotsOpened[0]) THEN RETURN, 1 ; nothing.. (press already button set)
	IF (opt_plotsOpened EQ 2) AND (id EQ wbttn_plotsOpened[1]) THEN RETURN, 2 ; nothing.. (press already button set)

	IF (opt_plotsOpened NE 1) AND (id EQ wbttn_plotsOpened[0]) THEN BEGIN
		WIDGET_CONTROL, wbttn_plotsOpened[0], SET_BUTTON=1
		WIDGET_CONTROL, wbttn_plotsOpened[1], SET_BUTTON=0
		opt_plotsOpened = 1
	ENDIF
	IF (opt_plotsOpened NE 2) AND (id EQ wbttn_plotsOpened[1]) THEN BEGIN
		WIDGET_CONTROL, wbttn_plotsOpened[0], SET_BUTTON=0
		WIDGET_CONTROL, wbttn_plotsOpened[1], SET_BUTTON=1
		opt_plotsOpened = 2
	ENDIF

	RETURN, opt_plotsOpened

END

;**************************************************************************************************
;**************************************************************************************************


FUNCTION fwdcemri_Parameters_write, ST_WIDGETS=st_widgets, ST_PARAM=st_param


	str_T10_tissue    = STRTRIM(STRING(FIX(st_param.T10_tissue*1000)),2)
	str_T10_blood     = STRTRIM(STRING(FIX(st_param.T10_blood*1000)),2)
	str_TR            = STRTRIM(STRING(FIX(st_param.tr*1000)),2)
	str_flip_angle    = STRTRIM(STRING(st_param.flip_angle_degrees,FORMAT='(F10.2)'),2)
	str_haematocrit   = STRTRIM(STRING(st_param.haematocrit,FORMAT='(F10.2)'),2)
	str_dose          = STRTRIM(STRING(st_param.dose, FORMAT='(F10.2)'),2)
	str_r1            = STRTRIM(STRING(st_param.r1 , FORMAT='(F10.2)'),2)
	str_frame_period    = STRTRIM(STRING(st_param.frame_period*60 , FORMAT='(F10.2)'),2) ; en segundos
	str_frame_injection = STRTRIM(STRING(FIX(st_param.nframe_injection)+1),2) ; Se suma 1, indice a 1 solo en el text
	str_nframes         = STRTRIM(STRING(FIX(st_param.nframes)),2)
	str_nframes_auc     = STRTRIM(STRING(FIX(st_param.nframes_auc)),2)
	str_nframes_adjust  = STRTRIM(STRING(FIX(st_param.nframes_adjust)),2)

	str_ktrans_rr      = STRTRIM(STRING(st_param.ktrans_rr , FORMAT='(F10.4)'),2)
	str_kep_rr         = STRTRIM(STRING(st_param.kep_rr ,     FORMAT='(F10.4)'),2)


	WIDGET_CONTROL, st_widgets.wtext_frame_period,      SET_VALUE=str_frame_period
	WIDGET_CONTROL, st_widgets.wtext_frame_injection,   SET_VALUE=str_frame_injection
	WIDGET_CONTROL, st_widgets.wtext_nframes,     SET_VALUE = str_nframes
	WIDGET_CONTROL, st_widgets.wtext_nframes_auc, SET_VALUE = str_nframes_auc
	;WIDGET_CONTROL, st_widgets.wtext_nframes_adjust, SET_VALUE = str_nframes_adjust
	WIDGET_CONTROL, st_widgets.wtext_dose , SET_VALUE=  str_dose
	WIDGET_CONTROL, st_widgets.wtext_r1,  SET_VALUE=  str_r1
	WIDGET_CONTROL, st_widgets.wtext_t10_tissue, SET_VALUE=  str_T10_tissue
	WIDGET_CONTROL, st_widgets.wtext_t10_blood, SET_VALUE=  str_T10_blood
	WIDGET_CONTROL, st_widgets.wtext_tr,  SET_VALUE=  str_TR

	WIDGET_CONTROL, st_widgets.wtext_haematocrit, SET_VALUE=  str_haematocrit
	WIDGET_CONTROL, st_widgets.wtext_flip_angle, SET_VALUE=  str_flip_angle

	WIDGET_CONTROL, st_widgets.wtext_ktrans_rr, SET_VALUE=  str_ktrans_rr
	WIDGET_CONTROL, st_widgets.wtext_kep_rr, SET_VALUE=  str_kep_rr


	RETURN, 1
END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_Parameters_read, ST_WIDGETS=st_widgets, ST_PARAM=st_param, ERROR_ID=error_id, ERROR_STR=error_str


;; PENDIENTE: limitar los valores válidos (negativos, valores no numéricos, etc)
	error_id = 0

	WIDGET_CONTROL, st_widgets.wtext_frame_period,    	GET_VALUE=str_frame_period
	WIDGET_CONTROL, st_widgets.wtext_frame_injection, 	GET_VALUE=str_frame_injection
	WIDGET_CONTROL, st_widgets.wtext_nframes,     		GET_VALUE=str_nframes
	WIDGET_CONTROL, st_widgets.wtext_nframes_auc, 		GET_VALUE=str_nframes_auc
	WIDGET_CONTROL, st_widgets.wtext_dose , 			GET_VALUE=str_dose
	WIDGET_CONTROL, st_widgets.wtext_r1, 				GET_VALUE=str_r1
	WIDGET_CONTROL, st_widgets.wtext_T10_tissue, 		GET_VALUE=str_T10_tissue
	WIDGET_CONTROL, st_widgets.wtext_T10_blood, 		GET_VALUE=str_T10_blood

	WIDGET_CONTROL, st_widgets.wtext_flip_angle, GET_VALUE=str_flip_angle
	WIDGET_CONTROL, st_widgets.wtext_tr, GET_VALUE=str_TR

	WIDGET_CONTROL, st_widgets.wtext_ktrans_rr, GET_VALUE=str_ktrans_rr
	WIDGET_CONTROL, st_widgets.wtext_kep_rr,    GET_VALUE=str_kep_rr

	IF str_t10_tissue NE 'T10 map' THEN BEGIN ;menos el caso de que haya un mapa de T10
		st_param.t10_tissue = FLOAT(str_T10_tissue[0])/1000.0  ; en segundos (está escrito en milisegundos)
	END
	st_param.t10_blood      = FLOAT(str_T10_blood[0])/1000.0   ; en segundos (está escrito en milisegundos en el campo)
	st_param.tr             = FLOAT(str_TR[0])/1000.0
	st_param.dose           = FLOAT(str_dose[0])
	st_param.r1             = FLOAT(str_r1[0])
	st_param.flip_angle_degrees = 0 > FLOAT(str_flip_angle[0])

	st_param.frame_period   = (0.001/60D) > FLOAT(str_frame_period[0])/60 ; en min-1

	st_param.nframes          = 1 > FIX(str_nframes[0])
	st_param.nframe_injection = 0 > (FIX(str_frame_injection[0])-1) < (st_param.nframes-1)
	st_param.nframes_adjust   = st_param.nframes - st_param.nframe_injection
	st_param.nframes_auc      = 1 > FIX(str_nframes_auc[0]) < (st_param.nframes_adjust-1)

	st_param.ktrans_rr      = FLOAT(str_ktrans_rr[0]) ; min-1
	st_param.kep_rr         = FLOAT(str_kep_rr[0])
	st_param.ve_rr          = st_param.ktrans_rr/st_param.kep_rr

	IF  st_param.nframe_injection LT 0 THEN BEGIN
		error_str = 'Wrong value: Injection frame' & error_id = -1 &	RETURN, -1
	ENDIF
	IF  st_param.nframe_injection GE (st_param.nframes-1) THEN BEGIN
		error_str = 'Wrong value: Injection frame' & error_id = -1 &	RETURN, -1
	ENDIF
	IF  st_param.nframes_auc GT st_param.nframes_adjust THEN BEGIN
		error_str = 'Wrong value: Frames for IAUC' & error_id = -1 &	RETURN, -1
	ENDIF

	RETURN, st_param
END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_ParametricAIF_write, ST_WIDGETS=st_widgets, ST_AIF=st_aif

	str_m1            = STRTRIM(STRING(st_aif.m1 , FORMAT='(F10.4)'),2)
	str_m2            = STRTRIM(STRING(st_aif.m2 , FORMAT='(F10.4)'),2)
	str_a1            = STRTRIM(STRING(st_aif.a1 , FORMAT='(F10.4)'),2)
	str_a2            = STRTRIM(STRING(st_aif.a2 , FORMAT='(F10.4)'),2)
	WIDGET_CONTROL, st_widgets.wtext_m1,  SET_VALUE=  str_m1
	WIDGET_CONTROL, st_widgets.wtext_m2,  SET_VALUE=  str_m2
	WIDGET_CONTROL, st_widgets.wtext_a1,  SET_VALUE=  str_a1
	WIDGET_CONTROL, st_widgets.wtext_a2,  SET_VALUE=  str_a2
	RETURN, 1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_ParametricAIF_read, ST_WIDGETS=st_widgets, ST_AIF=st_aif, ERROR_ID=error_id, ERROR_STR=error_str

	WIDGET_CONTROL, st_widgets.wtext_m1, GET_VALUE=str_m1
	WIDGET_CONTROL, st_widgets.wtext_m2, GET_VALUE=str_m2
	WIDGET_CONTROL, st_widgets.wtext_a1, GET_VALUE=str_a1
	WIDGET_CONTROL, st_widgets.wtext_a2, GET_VALUE=str_a2

	st_aif.model = 'biexponential'
	st_aif.m1  = FLOAT(str_m1[0])   ;  en min-1
	st_aif.m2  = FLOAT(str_m2[0])   ;  en min-1
	st_aif.a1  = FLOAT(str_a1[0])
	st_aif.a2  = FLOAT(str_a2[0])

	RETURN, st_aif

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_PlotROI_set, Pst, DATA=data

	(*Pst).vi1.Opolyline_ROI ->SetProperty, DATA=data, HIDE=0
	WIDGET_CONTROL, (*Pst).wd.wDraw1,  GET_VALUE = window_draw
	window_draw -> Draw, (*Pst).vi1.oView

	(*Pst).vi2.Opolyline_ROI ->SetProperty, DATA=data, HIDE=0
	WIDGET_CONTROL, (*Pst).wd.wDraw2,  GET_VALUE = window_draw
	window_draw -> Draw, (*Pst).vi2.oView

END

;**************************************************************************************************
;**************************************************************************************************
FUNCTION fwdcemri_Print_Data_3D, id_text, IMAG3D=imag3d

    total_value = TOTAL(imag3d)
    max_value = MAX(imag3d, MIN= min_value)
    tam = SIZE(imag3d, /DIMENSIONS)

    WIDGET_CONTROL, id_text, GET_VALUE=strarr_text1

    strarr_text1[0] = 'Total   3D: ' + STRTRIM(STRING(total_value, FORMAT='(F20.2)'),2)
    strarr_text1[1] = 'Max-Min 3D: [' + STRTRIM(STRING(max_value,  FORMAT='(F20.2)'),2) + $
        ', ' + STRTRIM(STRING(min_value, FORMAT='(F20.2)'),2) + ']'
    strarr_text1[2] = 'Size: [' + STRTRIM(tam[0],2) + ', ' + STRTRIM(tam[1],2) + ', ' + $
       STRTRIM(tam[2],2) + ']'
    WIDGET_CONTROL, id_text, SET_VALUE=strarr_text1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_Print_Data_2D, id_text, IMAG2D=imag2d

    total_value2d = TOTAL(imag2d)
    max_value2d = MAX(imag2d,MIN=min_value2d)

    WIDGET_CONTROL, id_text, GET_VALUE=strarr_text1

    strarr_text1[0] = 'Total   2D: ' + STRTRIM(STRING(total_value2d, FORMAT='(F20.2)'),2)
    strarr_text1[1] = 'Max-Min 2D: [' + STRTRIM(STRING(max_value2d,    FORMAT='(F20.2)'),2) + $
        ', ' + STRTRIM(STRING(min_value2d, FORMAT='(F20.2)'),2) + ']'
    WIDGET_CONTROL, id_text, SET_VALUE=strarr_text1

    RETURN, 1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_Redraw, Pst

(*Pst).vi1.Opolyline_ROI -> SetProperty,    HIDE=(*Pst).flag_drawROI EQ 0
(*Pst).vi2.Opolyline_ROI -> SetProperty,    HIDE=(*Pst).flag_drawROI EQ 0

WIDGET_CONTROL, (*Pst).wd.wDraw1,  GET_VALUE = window_draw
window_draw -> Draw, (*Pst).vi1.oView

WIDGET_CONTROL, (*Pst).wd.wDraw2,  GET_VALUE = window_draw
window_draw -> Draw, (*Pst).vi2.oView


RETURN, 1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_DrawInSameSize, Pst

	IF WIDGET_INFO((*Pst).wd.wbases_images[1], /MAP) THEN BEGIN
    	sizeb1 = WIDGET_INFO((*Pst).wd.wbases_images[0], /GEOMETRY)
		WIDGET_CONTROL, (*Pst).wd.wbases_images[1], XSIZE=sizeb1.xsize, YSIZE=sizeb1.ysize
	ENDIF

RETURN, 1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_RemakeIndexWithoutOutliers, Pst, MASK_LR=mask_lr

; Corrección 12-2011, elimina outliers de los indices de alta resolución

	range    = REFORM((*Pst).results.scales[*,(*Pst).opt_current_param,0]) ; rango válido de la escala del parámetro actual
	identifier = (*(*Pst).modelinfo).ids[(*Pst).opt_current_param]       ; parámetro actual, servirá para seleccionar los puntos que se eliminan por paleta

	arr_val = (*(*Pst).results.dataptr[identifier,0])
	pos_valid   = WHERE(arr_val GE range[0] AND arr_val LE range[1], pvalid)
	pos_invalid = WHERE(arr_val LT range[0] OR arr_val GT range[1],  pinvalid)
	IF pvalid EQ 0 THEN BEGIN
		ok = DIALOG_MESSAGE('No points in the ROI with the selected scale range', /INFO)
		RETURN, -1
	ENDIF

	mask_lr = *(*Pst).roi.mask_lr; máscara de baja resolución
	pos_lr  = WHERE(mask_lr EQ 1,ctp)
	IF pinvalid GT 0 THEN BEGIN ; ¡borra por outliers!
		pos_delete = pos_lr[pos_invalid]
		mask_lr[pos_delete] = 0
		; calcula los puntos de alta resolución a partir de la máscara de baja
		mask_points = ImIdxs_From_Indexes(MASK_LR=mask_lr,DIMENSIONS=(*Pst).inf.size_im,$
			ROI_RESOLUTION=(*Pst).roi.np_lr, MIN_P = (*Pst).roi.pr_ini, IDX_POINTS=idx_points)
		RETURN ,idx_points
	ENDIF ELSE BEGIN
		RETURN, *(*Pst).roi.idx_points
	ENDELSE

RETURN, *(*Pst).roi.idx_points


END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_ROIkinetics, Pst, PLOT=plot, IDX_POINTS=idx_points, WIN_NUMBER=win_number, XPOS=xpos, YPOS=ypos, TYPE_PLOT=type_plot,$
			ABSOLUTE_VALUE=absolute_Value

	; Absolute value - Raw MRI data (to export AIF, for example), instead of relative values

	;im =  REFORM((*(*Pst).data)[*,*,*,*])

	opt_absolutevalue = KEYWORD_SET(absolute_value)

	IF N_ELEMENTS(type_plot) EQ 0 THEN type_plot = 1
	CASE type_plot OF
	1 : BEGIN
		no_plot_sym = 1 & plot_error = 0
	END
	2 : BEGIN
		no_plot_sym = 1 & plot_error = 1
	END
	3 : BEGIN
		no_plot_sym = 0 & plot_error = 1
	END
	ENDCASE

	arr_dyn  = fwdcemri_DynDataFromIdx(Pst, IDX_POINTS=idx_points)
	n_points = (SIZE(arr_dyn,/DIMENSIONS))[0]

	arr_mean      = TOTAL(arr_dyn,1)/n_points

	IF opt_absolutevalue EQ 0 THEN BEGIN
		mean_baseline = MEAN(arr_mean[0:(*Pst).par.nframe_injection])
		scale_factor = 100.0
		arr_mean = (arr_mean/mean_baseline)*scale_factor
	ENDIF ELSE BEGIN
		mean_baseline = 1d & scale_factor = 1d
	ENDELSE

	arr_stddev = FLTARR((*Pst).inf.nframes)
	IF n_points GT 1 THEN BEGIN
		FOR i=0l, (*Pst).inf.nframes-1 DO BEGIN
			arr_stddev[i] = STDDEV(arr_dyn[*,i]/mean_baseline*scale_factor)
		ENDFOR
	ENDIF ELSE BEGIN
		arr_stddev = FLTARR((*Pst).inf.nframes)
	ENDELSE

	data_kinetics = [[*(*Pst).results.arr_time],[arr_mean],[arr_stddev]]

	IF KEYWORD_SET(plot) THEN BEGIN
		ok = plot_profile_KineticData(data_kinetics, WIN_NUMBER=win_number, PLOT_INJECTION=1, XPOS=xpos, YPOS=ypos,$
			INJECTION_POINT=(*Pst).par.nframe_injection,SCALE_FACTOR=mean_baseline/scale_factor, $
			NO_PLOT_LINE=0, NO_PLOT_PSYM=no_plot_sym, PLOT_ERROR=plot_error, TYPE_DATA=[3], PLOT_BASELINE=1)
	ENDIF
	data_kinetics[*,0]*=60 ; time in seconds
	data_kinetics = ROTATE(data_kinetics,4)
	RETURN, data_kinetics

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_ROI_sizeresolution, pst

	WIDGET_CONTROL, (*Pst).wd.wtext_ROI_size[0], GET_VALUE=str_ROI_size_x
	WIDGET_CONTROL, (*Pst).wd.wtext_ROI_size[1], GET_VALUE=str_ROI_size_y
	opt_ROIresolution = WIDGET_INFO((*Pst).wd.wdropl_ROIresolution, /DROPLIST_SELECT)
	roi_resolution = FIX((*Pst).wd.wdropl_ROIresolution_Array[opt_ROIresolution])

	roi_size = INTARR(2)
	roi_size[0]  = 1 > ((FIX(str_ROI_size_x))[0]) < ((*Pst).inf.size_im[0])
	roi_size[1]  = 1 > ((FIX(str_ROI_size_y))[0]) < ((*Pst).inf.size_im[1])

	FOR i=0,1 DO BEGIN
		IF roi_size[i] LT roi_resolution THEN roi_size[i] = roi_resolution
		IF roi_size[i] MOD roi_resolution NE 0 THEN roi_size[i] = roi_size[i] + (roi_resolution-(roi_size[i] MOD roi_resolution))
		IF roi_size[i] GT (*Pst).inf.size_im[i] THEN roi_size[i]-=roi_resolution
		WIDGET_CONTROL, (*Pst).wd.wtext_ROI_size[i],       SET_VALUE=STRTRIM(STRING(roi_size[i]),2)
	ENDFOR
	;-----------------------------------------------------
	(*Pst).inf.size_roi  =  roi_size
	(*Pst).roi.np_lr     =  roi_resolution
	;-----------------------------------------------------
	RETURN, 1

END
;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_ResetResults, pst


	FOR i=0l, N_ELEMENTS((*Pst).results.dataptr)/2-1 DO BEGIN
		undefine, *(*Pst).results.dataptr[i,0]
		undefine, *(*Pst).results.dataptr[i,1]
	ENDFOR

	undefine, *(*Pst).results.RCE_im

	undefine, *(*Pst).results.arr_Cp
	undefine, *(*Pst).results.arr_Signal_AIF
	undefine, *(*Pst).results.arr_signal_RR
	undefine, *(*Pst).results.arr_time
	undefine, *(*Pst).results.data_CT
	undefine, *(*Pst).results.arr_CT_RR
	undefine, *(*Pst).results.data_Signal_lr

	RETURN, 1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_ResetROI, pst

	(*Pst).vi2.Opolyline_AUX -> SetProperty, HIDE=1, DATA=[-1,-1]
	(*Pst).vi1.Opolyline_AUX -> SetProperty, HIDE=1, DATA=[-1,-1]
	(*Pst).vi2.Opolyline_ROI -> SetProperty, HIDE=(*pst).flag_drawROI, DATA=[-1,-1]
	(*Pst).vi1.Opolyline_ROI -> SetProperty, HIDE=(*pst).flag_drawROI, DATA=[-1,-1]

	(*Pst).flag_freeROIclose      = 0
	(*Pst).flag_freeROIinit       = 0
	(*Pst).flag_freeROIcontinuous = 0

	RETURN, 1

END

;****************************************************************************
;****************************************************************************

FUNCTION fwdcemri_Reshape_Windows, Pst, TWO_WINDOWS=two_windows

; Reshape window when a new image is loaded (not with zoom) and recalculates colorbar position
;
; Images can be not squares, so the relation X-Y can change to put the same colorbar size.
;
; The reserved size for colorbar is 25% of Xsize of image

opt_two_windows = KEYWORD_SET(two_windows)

tam = SIZE(*(*Pst).data, /DIMENSIONS)
xsize = tam[0]*1d
ysize = tam[1]*1d

(*Pst).opt.draw_size   = [(*Pst).opt.size_view_ini[0], ROUND((*Pst).opt.size_view_ini[0]*(ysize/xsize))]
(*Pst).opt.draw_size[1]+=((*Pst).opt.size_view_ini[0]*(*Pst).opt.colorbar_relativesize[0])  ; sum the reserved size for colorbar, independent of ysize (only image xsize accounts)

(*Pst).opt.draw_scroll_size = (*Pst).opt.draw_size

; Now we limit the draw size to avoid very big or very small windows in the screen

(*Pst).opt.draw_size[0] =  (*Pst).opt.view_size_limits_x[0] > (*Pst).opt.draw_size[0] < (*Pst).opt.view_size_limits_x[1]
(*Pst).opt.draw_size[1] =  (*Pst).opt.view_size_limits_y[0] > (*Pst).opt.draw_size[1] < (*Pst).opt.view_size_limits_y[1]

;--------------------------------------------------------------------
value_1 = ((*Pst).opt.draw_scroll_size[0]*1.0d/(*Pst).opt.draw_scroll_size[1])
value_2 = (xsize*1d/ysize)
scale_y = value_1/value_2

;scale_y = ysize/(ysize + xsize*(*Pst).opt.colorbar_relativesize[0]) ; se puede comprobar, es lo mismo menos el redondeo de draw_size..
;--------------------------------------------------------------------
;---------------------------------------------------------------------
; Viewplane_rect, coordinates in pixels
xMin = 0 & yMin = 0
xMax = tam[0] & yMax = tam[1] ; Para pasar a mm, cambiar n_1 por size_mm
xSpan = xMax-xMin
ySpan = yMax-yMin
xSpanmax = tam[0]-xMin
ySpanmax = tam[1]-yMin
val_frame_x = (*Pst).opt.val_frame1
val_frame_y = val_frame_x*scale_y*tam[0]/tam[1]  ; to compensate the scale in Y (model), and have the same frame in X and Y

ViewPlane_Rect = [xMin-val_frame_x*xSpanmax, yMin-val_frame_y*ySpanmax, $
   	xSpanmax*(1+val_frame_x*2), ySpanmax*(1+val_frame_y*2)] ; Con espacio alrededor

;---------------------------------------------------------------------
val_y_begin = 1 - ((1-scale_y)/2.2)
val_y_end   = 1 - ((1-scale_y)/100.0)

(*Pst).opt.position_colorbar_ini[1] = val_y_begin
(*Pst).opt.position_colorbar_ini[3] = val_y_end

position_colorbar = FLTARR(4)
position_colorbar[0] = (*Pst).opt.position_colorbar_ini[0]*xsize
position_colorbar[2] = (*Pst).opt.position_colorbar_ini[2]*xsize
position_colorbar[1] = (*Pst).opt.position_colorbar_ini[1]*ysize
position_colorbar[3] = (*Pst).opt.position_colorbar_ini[3]*ysize
;---------------------------------------------------------------------


(*Pst).vi1.omodel->setproperty, HIDE=1
(*Pst).vi1.omodel->reset
(*Pst).vi1.omodel->scale, 1, scale_y, 1
(*Pst).vi1.oview -> setProperty, VIEWPLANE_RECT = ViewPlane_Rect, HIDE=0
(*Pst).vi1.ohcolorbar -> setproperty, Position=position_colorbar
WIDGET_CONTROL, (*Pst).wd.wdraw1, XSIZE=(*Pst).opt.draw_size[0],  YSIZE=(*Pst).opt.draw_size[1], DRAW_XSIZE=(*Pst).opt.draw_scroll_size[0], DRAW_YSIZE=(*Pst).opt.draw_scroll_size[1]

(*Pst).vi1.omodel->setproperty, HIDE=0

;---------------------------------------------------------------------
IF opt_two_windows THEN BEGIN
	(*Pst).vi2.omodel->setproperty, HIDE=1
	(*Pst).vi2.omodel->reset
	(*Pst).vi2.omodel->scale, 1, scale_y, 1
	(*Pst).vi2.oview -> setProperty, VIEWPLANE_RECT = ViewPlane_Rect, HIDE=0
	(*Pst).vi2.ohcolorbar -> setproperty, Position=position_colorbar
	WIDGET_CONTROL, (*Pst).wd.wdraw2, XSIZE=(*Pst).opt.draw_size[0],  YSIZE=(*Pst).opt.draw_size[1], DRAW_XSIZE=(*Pst).opt.draw_scroll_size[0], DRAW_YSIZE=(*Pst).opt.draw_scroll_size[1]
	(*Pst).vi2.omodel->setproperty, HIDE=0

ENDIF

;---------------------------------------------------------------------
RETURN, 1

END

;****************************************************************************
;****************************************************************************

FUNCTION fwdcemri_SentEvents_AdjustPalette, Pst

; dibuja las rois y actualiza las paletas de color

WIDGET_CONTROL, (*Pst).wd.wtext_palette_maxval, SEND_EVENT={ID:0L, TOP:0L, HANDLER:0L}
RETURN,1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemri_TextWrite_pointval, id_text, point, VALUE=value, UNITS=units


IF N_ELEMENTS(value) NE 0 THEN BEGIN
	IF SIZE(value,/TYPE) EQ 1 THEN $
		str_type_val = ' = ' + STRTRIM(FIX(value[0]),2) $
	ELSE $
		str_type_val = ' = ' + STRTRIM(value[0], 2)
ENDIF ELSE BEGIN
	str_type_val = ''
ENDELSE

n_pos = N_ELEMENTS(point)
str_pos = '['
FOR i=0l, n_pos-1 DO BEGIN
	str_pos+= STRTRIM(point[i],2)
	IF i NE (n_pos-1) THEN $
		str_pos+=', ' $
	ELSE $
		str_pos+=']'
ENDFOR

IF N_ELEMENTS(units) NE 0 THEN str_units = units ELSE str_units = ''

WIDGET_CONTROL, id_text, SET_VALUE= str_pos+str_type_val + str_units

RETURN, 1

END

;***************************************************************************************************
;***************************************************************************************************

FUNCTION fwdcemri_DrawingEvents, ev, pst

    WIDGET_CONTROL, (*Pst).wd.wDraw1,  GET_VALUE = window_draw_1
    WIDGET_CONTROL, (*Pst).wd.wDraw2,  GET_VALUE = window_draw_2

    possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL', 'EXPOSE' ]
    EventType = possibleEventTypes[ev.type]

	IF eventtype EQ 'SCROLL' THEN BEGIN
		WIDGET_CONTROL, (*Pst).wd.wDraw1, GET_DRAW_VIEW=draw_view
		PRINT, draw_view
		PRINT, ev.x, ev.y
	END

	uname = WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'WINDOW1' : BEGIN
		window_draw = window_draw_1
		error = window_draw -> PickData((*Pst).vi1.oview, (*Pst).vi1.omodel, [ev.x,ev.y], Coord_view)
		id_text = (*Pst).wd.wtext_info1
	END
	'WINDOW2' : BEGIN
		window_draw = window_draw_2
		error = window_draw -> PickData((*Pst).vi2.oview, (*Pst).vi2.omodel, [ev.x,ev.y], Coord_view)
		id_text = (*Pst).wd.wtext_info2
	END
	ENDCASE
	;------------------------------------------------------------------------------
    IF (coord_view[0] LT 0 ) OR (coord_view[1] LT 0) THEN BEGIN
    	WIDGET_CONTROL, id_text, SET_VALUE=''
    	RETURN,-1
    ENDIF
    IF (coord_view[0] GE (*Pst).inf.size_im[0]) OR (coord_view[1] GE (*Pst).inf.size_im[1]) THEN BEGIN
    	WIDGET_CONTROL, id_text, SET_VALUE=''
    	RETURN,-1
    ENDIF
    point_fix = FIX(coord_view[0:1])
	;------------------------------------------------------------------------------
	type_action = ''
	opt_fast = 0l
	;------------------------------------------------------------------------------

	CASE EventType OF
	'MOTION': BEGIN
		type_action = 'Write_position'
		IF ((*Pst).flag_drawROI EQ 1) AND ((*Pst).roi.opt_typeROI EQ 2) THEN BEGIN
			type_action = 'draw_ROI_free:continue'
			opt_fast = 0 ; rapido... CAMBIO FEBRERO 2012
			BREAK
		END
	END
	'DOWN' : BEGIN ; ev.press = 1, botón izquierdo, ev.press = 4, botón derecho
		type_action = 'Write_position_and_value'
		IF (ev.press EQ 1) AND ((*Pst).flag_drawROI EQ 1) AND ((*Pst).roi.opt_typeROI EQ 0) THEN $
			type_action = 'draw_ROI_box' ; type_action
		IF (ev.press EQ 1) AND ((*Pst).flag_drawROI EQ 1) AND ((*Pst).roi.opt_typeROI EQ 1) THEN $
			type_action = 'draw_ROI_full' ; type_action
		IF (ev.press EQ 1) AND ((*Pst).flag_drawROI EQ 1) AND ((*Pst).roi.opt_typeROI EQ 2) THEN $
			type_action = 'draw_ROI_free:seg' ;

		IF (ev.press EQ 4) AND ((*Pst).flag_drawROI EQ 1) AND ((*Pst).roi.opt_typeROI EQ 2) THEN $
			type_action = 'draw_ROI_free:close' ;
		IF (ev.press EQ 4) AND (((*Pst).flag_drawROI EQ 0) OR ((*Pst).roi.opt_typeROI NE 2)) THEN $
			type_action = 'draw_dynProfile_HighRes'
	END
	'UP' : BEGIN ; levanta la pulsación de ratón
		IF (ev.release EQ 1) AND ((*Pst).flag_drawROI EQ 1) AND ((*Pst).roi.opt_typeROI EQ 2) THEN BEGIN
			type_action = 'draw_ROI_free:stop'
		ENDIF
	END
	ELSE: RETURN,-1
	ENDCASE

	IF opt_fast EQ 0 THEN BEGIN
		WIDGET_CONTROL,(*Pst).wd.wslide_z,  GET_VALUE=n_z
    	WIDGET_CONTROL,(*Pst).wd.wslide_n,  GET_VALUE=n_n
		tam = [(*Pst).inf.size_im[0],(*Pst).inf.size_im[1],(*Pst).inf.nslices,(*Pst).inf.nframes]
	END

	CASE type_action OF
	;------------------------------------------------------------------------------
	'Write_position' : BEGIN
		point_v2 = get_point_YupConvention(point_fix, DIMENSIONS=(*Pst).inf.size_im)
		ok = fwdcemri_textWrite_pointval(id_text, point_v2)
	END
	;------------------------------------------------------------------------------
	'Write_position_and_value' : BEGIN
		CASE uname OF
		'WINDOW1' : BEGIN
			val = (*(*Pst).data)[point_fix[0], point_fix[1], n_z, n_n]
		END
		'WINDOW2' : BEGIN
			val = ((*(*Pst).results.RCE_im)[*,*,n_z])[point_fix[0], point_fix[1]]*100.0
			str_units = ' % (RCE)'
		END
		ENDCASE
		point_v2 = get_point_YupConvention(point_fix, DIMENSIONS=(*Pst).inf.size_im)
		ok = fwdcemri_textWrite_pointval(id_text, point_v2, VALUE=val, UNITS=str_units)
	END
	;------------------------------------------------------------------------------
	'draw_ROI_box' : BEGIN

		ok = fwdcemri_ROI_sizeresolution(Pst)

		point_fix[0] = 0 > (point_fix[0]) < ((*Pst).inf.size_im[0]-(*Pst).inf.size_roi[0])
		point_fix[1] = 0 > (point_fix[1]-(*Pst).inf.size_roi[1]+1) < ((*Pst).inf.size_im[1]-1)
		(*Pst).inf.p_ini = point_fix
		PRINT, point_fix
		ok = fwdcemri_boxROI_set(Pst, POINT_INI=point_fix, DRAWTWO=1, ID_TEXT=[(*Pst).wd.wtext_info1,(*Pst).wd.wtext_info2])
		*(*Pst).roi.idx_points = IdxPointsFromBox(POINT_INI=(*Pst).inf.p_ini, SIZE_ROI=(*Pst).inf.size_roi, DIMENSIONS=(*Pst).inf.size_im)
		(*Pst).flag_noresetroi = 0
		ok = fwdcemri_Redraw(Pst)

	END
	;------------------------------------------------------------------------------
	'draw_ROI_free:seg' : BEGIN

		IF  ((*Pst).flag_freeROIclose EQ 0) THEN BEGIN
			(*Pst).vi1.Opolyline_ROI ->GetProperty, DATA=dat
			IF dat[0] NE -1 THEN BEGIN  ; No es el primer punto
				dat = [[dat], [coord_view[0:1]]]
			ENDIF ELSE BEGIN ; inicializa la polilínea por primera vez
				dat = [[coord_view[0:1]], [coord_view[0:1]]]
			ENDELSE
		ENDIF ELSE BEGIN ; inicializa la polilínea después de un close anterior
			(*Pst).flag_freeROIclose = 0
			dat =[[coord_view[0:1]], [coord_view[0:1]]]
		ENDELSE
		(*Pst).flag_freeROIinit  = 1
		(*Pst).flag_freeROIcontinuous = 1
		(*Pst).vi1.Opolyline_AUX ->SetProperty, DATA=[-1,-1], HIDE=1
		(*Pst).vi2.Opolyline_AUX ->SetProperty, DATA=[-1,-1], HIDE=1
		ok = fwdcemri_plotROI_set(Pst, DATA=dat)
		(*Pst).flag_noresetroi = 0
		ok = fwdcemri_activate_bases(Pst, /FORBID_ANALYSIS)
	END
	;------------------------------------------------------------------------------
	'draw_ROI_free:continue' : BEGIN

		;point_v2 = get_point_YupConvention(point_fix, DIMENSIONS=(*Pst).inf.size_im)
		;ok = fwdcemri_textWrite_pointval(id_text, point_v2)

		CASE uname OF
		'WINDOW1' : BEGIN
			val = (*(*Pst).data)[point_fix[0], point_fix[1], n_z, n_n]
		END
		'WINDOW2' : BEGIN
			val = ((*(*Pst).results.RCE_im)[*,*,n_z])[point_fix[0], point_fix[1]]*100.0
			str_units = ' % (RCE)'
		END
		ENDCASE
		point_v2 = get_point_YupConvention(point_fix, DIMENSIONS=(*Pst).inf.size_im)
		ok = fwdcemri_textWrite_pointval(id_text, point_v2, VALUE=val, UNITS=str_units)

		IF (*Pst).flag_freeROIinit NE 1 THEN RETURN,-1
		IF (*Pst).flag_freeROIcontinuous EQ 1 THEN BEGIN
			!EXCEPT=0
			(*Pst).vi1.Opolyline_ROI ->GetProperty, DATA=data
			n_p = (SIZE(data, /DIMENSIONS))[1]
			distance = ABS(coord_view[0]-data[0,n_p-1]) + ABS(coord_view[1]-data[1,n_p-1])
			IF distance GE 0 THEN BEGIN
				data = [[data], [coord_view[0:1]]]
				(*Pst).vi1.Opolyline_ROI ->SetProperty, DATA=data
				window_draw_1 ->  Draw, (*Pst).vi1.oView
				(*Pst).vi2.Opolyline_ROI ->SetProperty, DATA=data
				window_draw_2 -> Draw, (*Pst).vi2.oView
				ok = fwdcemri_activate_bases(Pst, /OPEN_ROI)
			ENDIF
			!EXCEPT=1
		ENDIF ELSE BEGIN
			(*Pst).vi1.Opolyline_ROI ->GetProperty, DATA=data
			n_p = (SIZE(data, /DIMENSIONS))[1]
			dat_temp = [[data[*,n_p-1]], [coord_view[0:1]]]
			(*Pst).vi1.Opolyline_AUX ->SetProperty, DATA=dat_temp, HIDE=0
			(*Pst).vi2.Opolyline_AUX ->SetProperty, DATA=dat_temp, HIDE=0
			ok = fwdcemri_Redraw(Pst)
		ENDELSE
		(*Pst).flag_noresetroi = 0
		RETURN,-1
	END
	;------------------------------------------------------------------------------
	'draw_ROI_free:stop' : BEGIN
		(*Pst).flag_freeROIcontinuous = 0

	END
	;------------------------------------------------------------------------------
	'draw_ROI_free:close' : BEGIN

		(*Pst).vi1.Opolyline_ROI -> GetProperty, DATA=data

		IF (SIZE(data,/DIMENSIONS))[1] GE 3 AND (*Pst).flag_freeROIclose EQ 0 THEN BEGIN

			IF MIN(data[*,0] EQ data [*,1]) EQ 1 THEN BEGIN
				data = data[*,1:(SIZE(data,/DIMENSIONS))[1]-1]
			ENDIF
			infpuntX = [[data[0,*]],[data[0,0]]]   ;cierra el polígono
			infpuntY = [[data[1,*]],[data[1,0]]]
			data = [infpuntX, infpuntY]
			(*Pst).vi1.Opolyline_AUX ->SetProperty, DATA=[-1,-1], HIDE=1
			(*Pst).vi2.Opolyline_AUX ->SetProperty, DATA=[-1,-1], HIDE=1

			IF (SIZE(data,/DIMENSIONS))[1] LE 3 THEN BEGIN
				ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main, $
					'Can not close single line, please make volumetric ROI', /INFO)
				RETURN,-1
			ENDIF
			;----------------------------------------------------------------------------------------
			*(*Pst).roi.idx_points = POLYFILLV(REFORM(ROUND(data[0,*])),  REFORM(ROUND(data[1,*])), (*Pst).inf.size_im[0], (*Pst).inf.size_im[1])
			;----------------------------------------------------------------------------------------
			IF (*(*Pst).roi.idx_points)[0] EQ -1 THEN BEGIN
				ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main, 'No points in ROI', /INFO) & RETURN,-1
			ENDIF
			(*Pst).flag_freeROIclose      = 1
			(*Pst).flag_freeROIinit       = 0
			(*Pst).flag_freeROIcontinuous = 0
			ok = fwdcemri_plotROI_set(Pst, DATA=data)
			ok = fwdcemri_activate_bases(Pst, /CLOSE_ROI)
			; Cerrado aquí para compatibilidad con ROIs guardados en disco
			;----------------------------------------------------------------------------------------
			(*Pst).flag_noresetroi = 0
			ok = fwdcemri_Activate_bases(Pst, /ALLOW_ANALYSIS)
		ENDIF
	END
	;------------------------------------------------------------------------------
    'draw_dynProfile_HighRes' : BEGIN
    	arr_dyn = REFORM((*(*Pst).data)[point_fix[0],point_fix[1],n_z,*])

    	opt_drawprofile = 1
    	opt_gaussian    = 0
    	;(*Pst).par.nframe_injection  = 3l
    	; ultima frame que se utiliza para calcular la concentracion sin contraste
    	IF opt_gaussian THEN $
			arr_dyn = Filter_Gaussian_1D(arr_dyn, SIGMA=3, /FLOATING)

    	win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)

    	CASE opt_drawprofile OF
    	1 : BEGIN ; direct profile
			st_par = (*Pst).par
			ok= fwdcemri_Parameters_read (ST_WIDGETS=(*Pst).wdp, ST_PARAM=st_par, ERROR_ID=error_id, ERROR_STR=error_str)
			ok= fwdcemri_Parameters_write(ST_WIDGETS=(*Pst).wdp, ST_PARAM=st_par)
			(*Pst).par = st_par
			*(*Pst).results.arr_time = INDGEN((*Pst).par.nframes)*(*Pst).par.frame_period ; actualizes time


			ok = plot_profile_dynMRI(arr_dyn, ARR_TIME=*(*Pst).results.arr_time, WIN_NUMBER=win_number, HR_POINT=point_fix[0:1],$
					FRAME_INJECTION = (*Pst).par.nframe_injection, XPOS=xpos, YPOS=ypos)
		END
		2 : BEGIN ; profile of relative increment of susceptibility (DSC images) R2 = -ln(S(t)/S(0))
			;arr_dyn[0:100] = S_o
			S_o = MEAN(arr_dyn[0:(*Pst).par.nframe_injection])
			ct_rel = -ALOG(arr_dyn/S_o)
			ok = plot_profile(ct_rel, WIN_NUMBER=win_number, HR_POINT=point_fix[0:1], NO_PLOT_PSYM=0, XPOS=xpos, YPOS=ypos)
		END
		ELSE:
		ENDCASE
		opt_prueba = 0
		IF opt_prueba EQ 1 THEN BEGIN
		;----------------------------------------
			im = REFORM((*(*Pst).data)[*,*,n_z,*])
			tot_begin = TOTAL(im[*,*,0:20-1],3)/20.0
			tot_max = im[*,*,0]*0
			FOR i=0l,(*Pst).inf.nframes-1 DO BEGIN
				tot_max = tot_max > im[*,*,i]
			ENDFOR
			im_diff = tot_max/tot_begin

			viewg, im_diff GT 1.5
		ENDIF

    END
    ;------------------------------------------------------------------------------
    ELSE:
    ENDCASE

END

;**********************************************************************************************************************
;**********************************************************************************************************************

FUNCTION fwdcemri_Perform_ModelAnalysis, ST_RESULTS=st_results, DATA=data, T1MAP=t1map, ST_AIF=st_aif, ST_MODELINFO=st_modelinfo, $
	ST_PAR=st_par, ST_ROI=st_roi, TYPE_AIF=type_aif, ERROR_ID=error_id


	WIDGET_CONTROL, /HOURGLASS

	error_id = 1l
	n_frames       = (SIZE(data,/DIMENSIONS))[3]
	n_points_lr    = N_ELEMENTS(*st_roi.idx_points_lr)

	;---------------------------------------------------------------------------------------------------
	data_Signal_lr = FLTARR(n_points_lr, n_frames)
	FOR i=0l, n_frames-1 DO BEGIN
		data_Signal_lr[*,i] = fwdcemri_DataFromSubidx(ST_ROI=st_roi, REFORM((data)[*,*,st_roi.slicez,i]))
	ENDFOR
	;---------------------------------------------------------------------------------------------------

	;IF st_modelinfo.using_AIF THEN BEGIN
	;IF st_modelinfo.using_RR  THEN BEGIN

	IF st_modelinfo.using_t10 THEN BEGIN
		IF N_ELEMENTS(t1map) LT 2 THEN BEGIN
			str_T1_tissue_ms = STRTRIM(STRING(FIX(st_par.T10_tissue*1000)),2)
			message_additional = 'Using T10 uniform value for tissue (' +  str_T1_tissue_ms + ' ms )'
			t1map_lr = REPLICATE(st_par.T10_tissue, n_points_lr)
			;option_common_t10 = 1
		ENDIF ELSE BEGIN
			t1map_lr =	fwdcemri_DataFromSubidx(ST_ROI=st_roi, REFORM((t1map)[*,*,st_roi.slicez]))
			;option_common_t10 = 0
		ENDELSE
	ENDIF

	CASE st_modelinfo.typemodel OF
	1: BEGIN ; 'CURVE

		opt_mindt = 0
		IF opt_mindt EQ 1 THEN BEGIN

			;im_ok = Analyze_DSCMRI_array(data_Signal_lr, FRAME_INI=st_par.nframe_injection, FRAME_TIME=st_par.frame_period*60, /ADD_MEAN, $
			;	IM_CBV=params_cbv, IM_MTT=params_mtt, IM_CBF=params_cbf, CT_REL=ct_rel)
			;MIND-T, possible analisis DSC-MRI

			names_param = ['RCBF', 'RCBV', 'RMTT']
			im_param = [[params_cbf], [params_CBV], [params_mtt]]
			names_stat  = ['']

			*st_results.data_signal_rce = ct_rel

		ENDIF ELSE BEGIN
			;-------------------------------------------------------
			params_auc = Analyze_AUC_array(data_Signal_lr, FRAME_INI=st_par.nframe_injection, $
				FRAME_END=st_par.nframe_injection+st_par.nframes_auc, /RELATIVE, /ADD_MEAN)
			;-------------------------------------------------------
			params_ttm = Analyze_TTM_array(data_Signal_lr, FRAME_INI=st_par.nframe_injection, /ADD_MEAN)
			;-------------------------------------------------------
			RCE_slice = RelativeContrastEnhancement( REFORM(data[*,*,st_roi.slicez,*]), FRAME_INJECTION=st_par.nframe_injection, ABSOLUTEMINMAX=0)
			params_RCE = fwdcemri_DataFromSubidx(ST_ROI=st_roi, REFORM(RCE_slice))*100.0 ; en porcentaje
			params_RCE = [params_RCE, MEAN(params_RCE)]; mean value
			;; hacer que sea el RCE de la suma, no la media de los resultados de los píxeles
			;---------------------------------------------------------------------

			Data_signal_RCE = FLTARR(n_points_lr+1, n_frames)
			Data_signal_RCE[0:n_points_lr-1,*]  = RCE_signalRel_from_Signal(data_Signal_lr, FRAME_INJECTION=st_par.nframe_injection) ; only calculates RCE in % per one
			Data_signal_RCE[n_points_lr,*]      = RCE_signalRel_from_Signal(TOTAL(data_Signal_lr,1)/n_points_lr, FRAME_INJECTION=st_par.nframe_injection) ; mean value

			*st_results.data_signal_rce = Data_signal_RCE

			names_param = ['RCE', 'IAUC', 'TTM']
			im_param = [[params_rce], [params_auc], [params_ttm]]
			names_stat  = ['']
		ENDELSE

	END
	2 : BEGIN ; 'HOFFMANN'
		opt_method_hoffmann = 1
		;------------------------------------------
		ok = Estimate_RCE_Array_HoffmannModel(data_Signal_lr,$
				ST_DATA  = st_par, $
				ARR_AIF  = data_aif_model,  $
				DATA_ERR = data_err, $
				DATA_SIGNAL_RCE = data_signal_rce, $
				IM_PARAM = im_param, $
				IM_PCERROR=im_pcerror,$
				IM_STATISTICS=im_statistics, $
				QUIET=1, $
				METHOD=opt_method_hoffmann,  $
				STAT_RANGE=st_par.nframes_stats,$
				NAMES_PARAM = names_param,$
				NAMES_STAT  = names_stat, $
				ERROR_ID=error_id, $
				ERROR_STR=error_str)
		;-------------------------------------------------------------------------------------
		IF error_id NE 0 THEN BEGIN & ok = DIALOG_MESSAGE(error_str, /INFO) &	RETURN,-1 &	ENDIF
		;-------------------------------------------------------------------------------------
		*st_results.data_signal_rce = data_signal_rce
	END
	;------------------------------------------------------------------------------------------
	3 : BEGIN
		IF st_modelinfo.model EQ 'TOFTS' THEN  opt_method_tofts =1 ELSE opt_method_tofts = 2 ; option 2 also estimates vp (extended Tofts)

		CASE type_aif OF
		0 : BEGIN & ok = DIALOG_MESSAGE('Error:AIF function not defined',/INFO, /CENTER) & RETURN, -1 & END
		1 : BEGIN
			IF SIZE(st_aif, /TNAME) EQ 'STRUCT' THEN BEGIN
				IF st_aif.model NE 'biexponential' THEN BEGIN
					arr_Cp = ParametricCp_to_CurveCp(ST_AIF=st_aif, ARR_TIME=*st_results.arr_time, DOSE=st_par.dose, INJECTION_FRAME=st_par.nframe_injection)
					option_convol = 1
				ENDIF ELSE BEGIN
					st_cp_model = st_aif
					option_convol = 0
				ENDELSE
			ENDIF ELSE BEGIN
				ok = DIALOG_MESSAGE('Error:AIF model not defined',/INFO, /CENTER) & RETURN, -1
			ENDELSE
		END
		ELSE : BEGIN
			IF N_ELEMENTS(*st_results.arr_Signal_AIF) NE st_par.nframes THEN BEGIN
				ok = DIALOG_MESSAGE('Error:AIF function in Tofts model is not compatible with data, diferent number of samples',/INFO, /CENTER)
				RETURN, -1
			ENDIF

			S_o = MEAN((*st_results.arr_Signal_AIF)[0:st_par.nframe_injection])
			;--------------------------------------------------------------------------------------
			arr_Cb = RCE_Signal_to_tissueConcentration_Barboriak(*st_results.arr_Signal_AIF, S_o=S_o, T10=st_par.t10_blood, TR=st_par.TR, R1=st_par.r1, ANGLE=st_par.flip_angle_degrees*!dpi/180)
			arr_Cp = arr_Cb/(1.0-st_par.haematocrit)
			;arr_Cp[0:st_par.nframe_injection] = 0 ; by definition
			;--------------------------------------------------------------------------------------
			option_convol=1
		END
		ENDCASE
		IF SIZE(message_additional, /TNAME)  EQ 'STRING' THEN $
			ok = DIALOG_MESSAGE(message_additional, /INFO, /CENTER)
		;------------------------------------------------------------------------------------------
		time_initial = SYSTIME(1)
		ok = Estimate_Ct_Array_ToftsModel(data_Signal_lr, t1map_lr,$
			ST_DATA=st_par, $
			ST_CP_MODEL=st_cp_model,$
			ARR_CP= arr_cp, $
			METHOD=opt_method_tofts, $
			IMAG_FLAG=im_flag, $
			IM_PARAM=im_param, $
			IM_PCERROR=im_pcerror,$
			IM_STATISTICS=im_statistics,$
			ARR_AUX=arr_aux, $
			DATA_ERR=data_err, $
			DATA_CT=data_CT,  $
			ARR_TIME=arr_t_tofts,   $
			STAT_RANGE=st_par.nframes_stats,$
			COMMON_T10=option_common_t10,$
			ERROR_ID=error_id,$
			ERROR_STR=error_str,$
			NAMES_PARAM = names_param,$
			NAMES_STAT  = names_stat, $
			ARR_STATUS  = arr_status,$
			OPTION_CONVOL=option_convol)
		;--------------------------------------------------------------------------------------------
		IF error_id NE 0 THEN BEGIN & ok = DIALOG_MESSAGE(error_str, /INFO) &	RETURN,-1 &	ENDIF
		PRINT, arr_status
		time_elapsed = Get_elapsedtime(time_initial)
		PRINT, 'Time elapsed:', time_elapsed
		;--------------------------------------------------------------------------------------------
		IF option_convol EQ 1 THEN BEGIN
			st_results.tofts_type = 1
			*st_results.arr_Cp = arr_cp
		ENDIF ELSE BEGIN
			st_results.tofts_type = 2
			;*st_results.arr_Cp = arr_cp  ; ¿?
			*st_results.cp_model = st_cp_model
		ENDELSE
		*st_results.data_ct     = data_ct
	END
	4 : BEGIN

		CASE type_aif OF
		0 : BEGIN & ok = DIALOG_MESSAGE('Error: AIF function not defined',/INFO, /CENTER) & RETURN, -1 & END
		1 :
		ELSE : BEGIN
			BEGIN & ok = DIALOG_MESSAGE('Error: AIF function must be defined analitically',/INFO, /CENTER) & RETURN, -1 & END
		END
		ENDCASE
		;--------------------------------------------------------------------------------------------
		ok = Estimate_RCE_array_LarssonModel(data_Signal_lr,$
			ST_DATA    =  st_par, $
			ST_CP_MODEL = st_cp_model, $
			DATA_ERR   =  data_err, $
			DATA_SIGNAL_RCE   = data_signal_rce, $
			IM_PARAM   = im_param, $
			IM_PCERROR= im_pcerror, $
			IM_STATISTICS=im_statistics, $
			QUIET=1, $
			STAT_RANGE=st_par.nframes_stats,$
			NAMES_PARAM=names_param,$
			NAMES_STAT =names_stat, $
			ERROR_ID=error_id,$
			ERROR_STR=error_str)
		;------------------------------------------------------------------------------------------
		IF error_id NE 0 THEN BEGIN & ok = DIALOG_MESSAGE(error_str, /INFO) &	RETURN,-1 &	ENDIF
		;-------------------------------------------------------------------------------------
		*st_results.data_signal_rce = data_signal_rce
	END
	5 : BEGIN
		;--------------------------------------------------------------------------------------------------------
		IF N_ELEMENTS(*st_results.arr_signal_RR) LE 1 THEN BEGIN
			ok = DIALOG_MESSAGE('Error: No reference region selected for RR model',	/INFO,	/CENTER)
			RETURN, -1
		ENDIF
		IF N_ELEMENTS(*st_results.arr_signal_RR) NE st_par.nframes THEN BEGIN
			ok = DIALOG_MESSAGE('Error:Reference region is not compatible with data, diferent number of samples',	/INFO,	/CENTER)
			RETURN, -1
		ENDIF
		IF SIZE(message_additional, /TNAME)  EQ 'STRING' THEN $
			ok = DIALOG_MESSAGE(message_additional, /INFO, /CENTER)
		;------------------------------------------------------------------------------------------
		ok = Estimate_Ct_Array_YankeelovModel(data_Signal_lr, t1map_lr,$
			ST_DATA=st_par, $
			ARR_SIGNAL_RR= *st_results.arr_signal_RR, $
			IMAG_FLAG=im_flag, $
			IM_PARAM=im_param, $
			IM_PCERROR=im_pcerror,$
			IM_STATISTICS=im_statistics,$
			ARR_AUX=arr_aux, $
			DATA_ERR=data_err, $
			ARR_CT_RR=arr_Ct_RR,   $
			DATA_CT_t=data_Ct,      $
			ARR_TIME=arr_t_tofts,   $
			STAT_RANGE=st_par.nframes_stats,$
			COMMON_T10=option_common_t10,$
			KTRANS_RR=st_par.ktrans_rr, KEP_RR=st_par.kep_rr,$
			ERROR_ID=error_id, $
			ERROR_STR=error_str,$
			NAMES_PARAM = names_param,$
			NAMES_STAT  = names_stat)
		;------------------------------------------------------------------------------------------
		IF error_id NE 0 THEN BEGIN & ok = DIALOG_MESSAGE(error_str, /INFO) &	RETURN,-1 &	ENDIF
		;------------------------------------------------------------------------------------------
		*st_results.arr_Ct_RR  = arr_Ct_RR
		*st_results.data_ct    = data_ct
	END
	ENDCASE
	im_param  = finite_values(im_param)
	IF N_ELEMENTS(im_pcerror) GT 0 THEN	$
		im_pcerror= finite_values(im_pcerror)
	;---------------------------------------------------------------------------------------------
	FOR i=0, N_ELEMENTS(st_results.dataptr)-1 DO BEGIN
		undefine, *st_results.dataptr[i]
	ENDFOR

	FOR i=0, st_modelinfo.nparam-1 DO BEGIN
		str_id = STRTRIM(st_modelinfo.ids[i],2)
		posp = WHERE(STRTRIM(names_param,2) EQ str_id[0], npp)
		poss = WHERE(names_stat  EQ str_id[0], nps)
		IF ((npp EQ 1) OR (nps EQ 1)) EQ 0 THEN BEGIN
			ok = DIALOG_MESSAGE('Internal error #paranob454, please contact software supplier',	/INFO, /CENTER)
			RETURN, -1
		ENDIF
		IF npp EQ 1 THEN BEGIN
			*st_results.dataptr[i,0]=im_param[*,posp]
			IF st_modelinfo.sd_valid[i] THEN $
				*st_results.dataptr[i,1]=im_pcerror[*,posp]
		ENDIF ELSE BEGIN
			*st_results.dataptr[i,0]=im_statistics[*,poss]
		ENDELSE
	ENDFOR
	;--------------------------------------------------------------------------------------------

	; "st_roi.pr_ini" es como el "st_inf.p_ini"  pero general para ROIS cuadradas, freees, etc...
	;  El resultado se hereda de st_inf.p_ini en el caso de ROIS cuadradas, y es independiente para otros casos
	;  - Lo mismo para pr_ini, roi_sz, y np_lr

	;------------------------------------------------------------------------------------------------
	*st_results.data_Signal_lr   = data_Signal_lr
	*st_results.data_Signal_ROI  = TOTAL(data_Signal_lr,1)/n_points_lr

	error_id = 0
	RETURN, st_results

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwBaseStruct_Window, PARENT_BASE=parent_base, MAIN_BASE=main_base, TITLE=title, ST_OPTIONS=st_options, SIZE_IMAGE=size_image, $
	ST_OBJECTS=st_objects , PALETTE_ROI=palette_roi, PALETTE_IMAGE=palette_image, DRAW_ID=draw_id, UNAME=uname, UVALUE=uvalue, ROI=roi

base_win    = WIDGET_BASE(parent_base, /COLUMN, FRAME=0, /ALIGN_CENTER)

IF N_ELEMENTS(title) EQ 0 THEN title=''
opt_roi = KEYWORD_SET(roi)

label_draw  = WIDGET_LABEL(base_win, VALUE=title)
base_draw   = WIDGET_BASE(base_win, /COLUMN, FRAME=0, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
draw        = WIDGET_DRAW(Base_draw, RETAIN=2, COLOR_MODEL=0, RENDERER=1,GRAPHICS_LEVEL=2, $
					XSIZE=st_options.draw_size[0], YSIZE=st_options.draw_size[1], $
					X_SCROLL_SIZE=st_options.draw_scroll_size[0], Y_SCROLL_SIZE=st_options.draw_scroll_size[1], SCROLL=1, $
					MOTION_EVENTS=1, BUTTON_EVENTS=1, TRACKING_EVENTS=0, VIEWPORT_EVENTS=1, EXPOSE_EVENTS=0,$
					UVALUE=uvalue, UNAME=uname)

WIDGET_CONTROL, draw, GET_VALUE  = ds_window


IF N_ELEMENTS(st_objects) NE 0 THEN BEGIN
	st_objects.opalette->GetProperty, RED=r1,GREEN=g1,BLUE=b1
ENDIF ELSE BEGIN
	IF N_ELEMENTS(palette_roi) GT 0 THEN LOADCT, palette_image
	TVLCT, r1,g1,b1, /GET
ENDELSE

IF N_ELEMENTS(size_image) EQ 0 THEN size_image = [100,100]


;--------------------------------------------
;ds_colorbar = OBJ_NEW( 'IDLgrColorbar', r,g,b, SHOW_AXIS=2, HIDE=0, DIMENSIONS = [256, 16])
ds_container   = OBJ_NEW('IDL_Container')
ds_scene       = OBJ_NEW('IDLgrScene')
ds_view        = OBJ_NEW('IDLgrView', COLOR=st_options.background_color_mri)
ds_model       = OBJ_NEW('IDLgrModel')
ds_palette     = OBJ_NEW('IDLgrPalette',r1,g1,b1)
ds_image       = OBJ_NEW('IDLgrImage', PALETTE=ds_palette,    GREYSCALE=0) ;INTERLEAVE=2 ;

ds_polyline_ROI  = OBJ_NEW('IDLgrPolyline', DATA=[-1,-1], HIDE=0, COLOR=[255,0,0], THICK=1)
ds_polyline_AUX  = OBJ_NEW('IDLgrPolyline', DATA=[-1,-1], HIDE=1, COLOR=[255,0,0], THICK=1, LINESTYLE=1)
ds_text          = OBJ_NEW('IDLgrText', 'File Size Exceeded')

ds_hcolorbar = OBJ_NEW('HColorBar2', Palette=ds_palette, $
	POSITION=st_options.position_colorbar_ini*size_image[0], $
	COLOR=st_options.colorbar_labelcolor, $
	FONTSIZE=st_options.colorbar_fontsize)
ds_text     -> SetProperty, HIDE=1

IF opt_roi THEN BEGIN

	IF N_ELEMENTS(st_objects) NE 0 THEN BEGIN
		st_objects.opaletteROI->GetProperty, RED=r2,GREEN=g2,BLUE=b2

		st_objects.opolyline_ROI->GetProperty, DATA=data
		ds_polyline_ROI->SetProperty, DATA=data

	ENDIF ELSE BEGIN
		IF N_ELEMENTS(palette_roi) GT 0 THEN LOADCT, palette_roi
		TVLCT, r2,g2,b2, /GET
	ENDELSE

	ds_paletteROI  = OBJ_NEW('IDLgrPalette', r2,g2,b2)
	ds_imageROI    = OBJ_NEW('IDLgrImage',   PALETTE=ds_paletteROI, GREYSCALE=0, HIDE=0, $
	ALPHA_CHANNEL=1.0, BLEND_FUNCTION=[3,4], INTERLEAVE=2) ;INTERLEAVE=2 ;
	; nueva imagen para la mascara auxiliar, con transparencia
	ds_imageMASK   = OBJ_NEW('IDLgrImage', GREYSCALE=0, HIDE=1, $
	ALPHA_CHANNEL=1.0, BLEND_FUNCTION=[3,4], INTERLEAVE=2)
ENDIF ELSE BEGIN

	ds_imageROI   = -1L
	ds_imageMASK  = -1L
	ds_paletteROI = -1L
ENDELSE

std = { $
    ocontainer : ds_container,$
    owindow : ds_window,$
    oscene  : ds_scene, $
    oview   : ds_view,  $
    omodel  : ds_model, $
    opalette: ds_palette,$
    opaletteROI: ds_paletteROI,$
    oimage     : ds_image,     $
    oimageROI  : ds_imageROI,  $
    oimageMASK : ds_imageMASK, $
    otext   : ds_text,$
    opolyline_ROI : ds_polyline_ROI, $
    opolyline_AUX : ds_polyline_AUX, $
    ohcolorbar    : ds_hcolorbar     $
}

std.oscene    -> Add, std.oview
std.omodel    -> Add, std.oimage
std.omodel    -> Add, std.otext
std.omodel    -> Add, std.opolyline_ROI
std.omodel    -> Add, std.opolyline_AUX
std.oview     -> Add, std.omodel
std.oview     -> Add, std.ohcolorbar ; its another model

std.ocontainer-> Add, [std.owindow, std.oscene, std.omodel, std.oview, std.opalette, $
	std.otext, std.opolyline_ROI, std.opolyline_AUX, std.ohcolorbar, std.oimage]

IF opt_roi EQ 1 THEN BEGIN
	std.omodel    -> Add, std.oimageROI
	std.omodel    -> Add, std.oimageMASK
	std.ocontainer-> Add, [std.opaletteROI, std.oimageROI, std.oimageMASK]
ENDIF

draw_id   = TEMPORARY(draw)
main_base = TEMPORARY(base_win)


RETURN, std


END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwBaseStruct_Frames, PARENT_BASE=parent_base, MAIN_BASE=main_base, FRAME=frame, SET_VALUES=set_values, ACTIVE=active


;-------------------------------------------------------------------
IF N_ELEMENTS(set_values) NE 0 THEN BEGIN
	current_slide_n = set_values[0]
	current_slide_z = set_values[1]
	max_slide_n     = set_values[2]-1
	max_slide_z     = set_values[3]-1
	value_label_z = STRTRIM(STRING(current_slide_z),2)
	value_label_n = STRTRIM(STRING(current_slide_n),2)
ENDIF ELSE BEGIN
	current_slide_z = 0 & max_slide_z = 1 & current_slide_n = 0 & max_slide_n = 1
	value_label_z = '' & value_label_n = ''
ENDELSE

IF N_ELEMENTS(active) EQ 0 THEN BEGIN
	opt_sensitive_n = 1 & opt_sensitive_z = 1
ENDIF ELSE BEGIN
	opt_sensitive_n = active[0] & opt_sensitive_z = active[1]
ENDELSE
;-------------------------------------------------------------------

base_1=WIDGET_BASE(parent_base, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER, YPAD=0, XPAD=0, FRAME=frame)
	base_1a=WIDGET_BASE(base_1, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER, YPAD=0)
	base_1b=WIDGET_BASE(base_1, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER, YPAD=0)

slide_n = WIDGET_SLIDER(base_1a, MAXIMUM=max_slide_n, MINIMUM=0, SCR_XSIZE=size_slider_1, $
    VALUE=current_slide_n, UVALUE='S_CHANGE', UNAME = 'SLICE_N',  SENSITIVE=opt_sensitive_n, $
    TITLE='Dynamic frame', /SUPPRESS_VALUE)
label_n = WIDGET_LABEL(base_1a, VALUE=value_label_n, SCR_XSIZE=25, SCR_YSIZE=25)

IF max_slide_z EQ 0 THEN max_slide_z++

slide_z = WIDGET_SLIDER(base_1b, MAXIMUM=max_slide_z, MINIMUM=0, SCR_XSIZE=size_slider_2, $
    VALUE=current_slide_z, UVALUE='S_CHANGE', UNAME = 'SLICE_Z',  SENSITIVE=opt_sensitive_z, $
    TITLE='Slice', /SUPPRESS_VALUE)
label_z= WIDGET_LABEL(base_1b, VALUE=value_label_z, SCR_XSIZE=20, SCR_YSIZE=20)

main_base = TEMPORARY(base_1)

st_widgets_frames = {slide_n : slide_n,$
					slide_z : slide_z,$
					label_n  : label_n,$
					label_z : label_z}

RETURN, st_widgets_frames

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwBaseStruct_Scales, PARENT_BASE=parent_base, MAIN_BASE=main_base, FRAME=frame

;---------------------------------------------------------------------

main_base = WIDGET_BASE(parent_base, /COLUMN, /ALIGN_LEFT, /BASE_ALIGN_LEFT, YPAD=0, XPAD=0, FRAME=frame)

base_1=WIDGET_BASE(main_base, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT, YPAD=0, XPAD=0)
		base_1a=WIDGET_BASE(base_1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT, YPAD=0)
		base_1b=WIDGET_BASE(base_1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT, YPAD=0)
base_2=WIDGET_BASE(main_base, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT, YPAD=0, XPAD=0)
       	base_2a=WIDGET_BASE(base_2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT, YPAD=0)


slide_minval = WIDGET_SLIDER(base_1a, MAXIMUM=200, MINIMUM=0, SCR_XSIZE=size_slider_3, $
    VALUE=100, UVALUE='VAL_CHANGE', UNAME = 'SLIDE_MINVAL',  $
    TITLE='min value (RCE %) ', /SUPPRESS_VALUE)
;label_minval   = WIDGET_TEXT(basse_a1, VALUE='', UVALUE='TEXT_MINVAL',/EDITABLE, XSIZE=3, SCR_YSIZE=18)
;label_minval = WIDGET_LABEL(base_1a, VALUE='', SCR_XSIZE=25, SCR_YSIZE=20)
text_minval  = WIDGET_TEXT(base_1a,  VALUE= '',  UVALUE='ADJUST_VALUES', $
	/EDITABLE, XSIZE=8, UNAME= 'VALUE_MINVAL')

slide_maxval = WIDGET_SLIDER(base_1b, MAXIMUM=1000, MINIMUM=0, SCR_XSIZE=size_slider_3, $
    VALUE=300, UVALUE='VAL_CHANGE', UNAME = 'SLIDE_MAXVAL',  $
    TITLE='max value (RCE %)', /SUPPRESS_VALUE)
;label_maxval  = WIDGET_LABEL(base_1b, VALUE='', SCR_XSIZE=25, SCR_YSIZE=20)
text_maxval  = WIDGET_TEXT(base_1b,  VALUE= '',  UVALUE='ADJUST_VALUES', $
	/EDITABLE, XSIZE=8, UNAME= 'VALUE_MAXVAL')

slide_threshold = WIDGET_SLIDER(base_2a, MAXIMUM=100, MINIMUM=0, SCR_XSIZE=size_slider_3, $
   VALUE=10, UVALUE='VAL_CHANGE', UNAME = 'SLIDE_THRESHOLD',  $
    TITLE='Threshold (min val %)', /SUPPRESS_VALUE)
label_threshold = WIDGET_LABEL(base_2a, VALUE='', SCR_XSIZE=25, SCR_YSIZE=20)


st_widgets_scales = {slide_minval: slide_minval,$
					slide_maxval : slide_maxval,$
					text_minval  : text_minval, $
					text_maxval  : text_maxval,  $
					slide_threshold:slide_threshold,$
					label_threshold :label_threshold }


RETURN, st_widgets_scales

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwBaseStruct_Zoom, PARENT_BASE=parent_base, MAIN_BASE=main_base, FRAME=frame


main_base=WIDGET_BASE(parent_base, /ROW,  /ALIGN_CENTER, /BASE_ALIGN_CENTER, FRAME=frame)
	base_1=WIDGET_BASE(main_base, /ROW,  /ALIGN_CENTER, /BASE_ALIGN_CENTER)
	base_2=WIDGET_BASE(main_base, /ROW,  /ALIGN_CENTER, /BASE_ALIGN_CENTER, /NONEXCLUSIVE)
bttn_zoommenos  = WIDGET_BUTTON(base_1, VALUE='Zoom -',   UNAME='ZOOM_MENOS', UVALUE='ZOOM_CHANGE')
bttn_zoommas    = WIDGET_BUTTON(base_1, VALUE='Zoom +',   UNAME='ZOOM_MAS',   UVALUE='ZOOM_CHANGE')
bttn_interp     = WIDGET_BUTTON(base_2, VALUE= 'Interp.', UVALUE='INTERPOLATION')

st_widgets_zoom = {base_main:main_base, $
			bttn_zoommenos:bttn_zoommenos, $
			bttn_zoommas:bttn_zoommas, $
			bttn_interp:bttn_interp}

RETURN, st_widgets_zoom

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwBaseStruct_Parameters, PARENT_BASE=parent_base, ST_PARAM=st_param, ST_AIF=st_aif, $
	ARR_INI=arr_ini, TYPE=type, UVALUE=uvalue

IF N_ELEMENTS(uvalue) EQ 0 THEN uvalue= 'CHOOSE_MODE'

base_x = WIDGET_BASE(parent_base, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
base_y = WIDGET_BASE(parent_base, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
	base_y1 = WIDGET_BASE(base_y, /COLUMN, /ALIGN_CENTER, /BASE_ALIGN_CENTER)

base_tab = WIDGET_TAB(base_y1, /ALIGN_CENTER, LOCATION=0)
	base_t1 = WIDGET_BASE(base_tab, COLUMN=2, /ALIGN_LEFT, /BASE_ALIGN_LEFT, TITLE='MR signal')
		base_t1a = WIDGET_BASE(base_t1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t1b = WIDGET_BASE(base_t1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t1c = WIDGET_BASE(base_t1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t1d = WIDGET_BASE(base_t1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t1e = WIDGET_BASE(base_t1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t1f = WIDGET_BASE(base_t1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
	base_t2 = WIDGET_BASE(base_tab, COLUMN=2, /ALIGN_LEFT, /BASE_ALIGN_LEFT, TITLE='Contrast Agent')
		base_t2a = WIDGET_BASE(base_t2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t2b = WIDGET_BASE(base_t2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t2c = WIDGET_BASE(base_t2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t2d = WIDGET_BASE(base_t2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t2e = WIDGET_BASE(base_t2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
	base_t3_1 = WIDGET_BASE(base_tab, COLUMN=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT, TITLE='AIF')

		base_t3_2 = WIDGET_BASE(base_t3_1, COLUMN=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t3   = WIDGET_BASE(base_t3_1, COLUMN=2, /ALIGN_LEFT, /BASE_ALIGN_LEFT)

		base_t3a = WIDGET_BASE(base_t3, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t3b = WIDGET_BASE(base_t3, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t3c = WIDGET_BASE(base_t3, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t3d = WIDGET_BASE(base_t3, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
	base_t4 = WIDGET_BASE(base_tab, COLUMN=2, /ALIGN_LEFT, /BASE_ALIGN_LEFT, TITLE='RR')
		base_t4a = WIDGET_BASE(base_t4, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t4b = WIDGET_BASE(base_t4, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t4c = WIDGET_BASE(base_t4, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t4d = WIDGET_BASE(base_t4, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)


str_ini = ''

text_t10_tissue  =  WIDGET_TEXT(base_t1a, VALUE=str_ini, XSIZE=6, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_t10_blood = WIDGET_LABEL(base_t1a, VALUE=' T10 in tissue (ms)')

text_t10_blood   =  WIDGET_TEXT(base_t1b, VALUE=str_ini, XSIZE=6, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_t10_blood = WIDGET_LABEL(base_t1b, VALUE=' T10 in blood (ms)')

text_tr  =  WIDGET_TEXT(base_t1c, VALUE=str_ini, XSIZE=6, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_tr = 	WIDGET_LABEL(base_t1c, VALUE=' TR (ms)')

text_frame_period  =  WIDGET_TEXT(base_t1d, VALUE=str_ini, XSIZE=6, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_frame_period = 	WIDGET_LABEL(base_t1d, VALUE=' Frame period (s) ')

text_nframes  =  WIDGET_TEXT(base_t1e, VALUE=str_ini, XSIZE=6,  $
	/ALIGN_LEFT, SENSITIVE=1, EDITABLE=0, UVALUE=uvalue)
label_nframes =  WIDGET_LABEL(base_t1e, VALUE=' Nº of frames')

text_flip_angle = WIDGET_TEXT(base_t1f, VALUE=str_ini, XSIZE=6, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE,  UVALUE=uvalue)
label_flip_angle = WIDGET_LABEL(base_t1f, VALUE=' Flip angle (Degrees)')


text_dose  =  WIDGET_TEXT(base_t2a, VALUE=str_ini, XSIZE=6, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_dose =  WIDGET_LABEL(base_t2a, VALUE=' Dose (mmol/kg)')

text_r1  =  WIDGET_TEXT(base_t2b, VALUE=str_ini, XSIZE=6, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_r1 = 	WIDGET_LABEL(base_t2b, VALUE=' Relaxivity (mM s-1) ')

text_haematocrit = WIDGET_TEXT(base_t2c, VALUE=str_ini, XSIZE=6, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_haematocrit = WIDGET_LABEL(base_t2c, VALUE=' Haematocrit ')

text_frame_injection  =  WIDGET_TEXT(base_t2d, VALUE=str_ini, XSIZE=4,  $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_frame_ini =  WIDGET_LABEL(base_t2d, VALUE=' Injection frame [1,n)')


text_nframes_auc  =  WIDGET_TEXT(base_t2e, VALUE=str_ini, XSIZE=4,  $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_nframes_auc =  WIDGET_LABEL(base_t2e, VALUE=' Frames for IAUC')

text_ktrans_rr  =  WIDGET_TEXT(base_t4a, VALUE=str_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_ktrans_rr  = WIDGET_LABEL(base_t4a, VALUE=' Ktrans(RR) (min-1)')

text_kep_rr   =  WIDGET_TEXT(base_t4b, VALUE=str_ini, XSIZE=8,$
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_kep_rr  = 	WIDGET_LABEL(base_t4b, VALUE='kep(RR) (min-1)')

label_parametric_values = WIDGET_LABEL(base_t3_2, VALUE=' Biexponential model: ')

text_m1  =  WIDGET_TEXT(base_t3a, VALUE=str_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_m1 = 	WIDGET_LABEL(base_t3a, VALUE=' m1 (min-1) ')

text_m2  =  WIDGET_TEXT(base_t3b, VALUE=str_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_m2 = 	WIDGET_LABEL(base_t3b, VALUE=' m2 (min-1) ')

text_a1  =  WIDGET_TEXT(base_t3c, VALUE=str_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_a1 = 	WIDGET_LABEL(base_t3c, VALUE=' a1 (kg/l) ')

text_a2  =  WIDGET_TEXT(base_t3d, VALUE=str_ini, XSIZE=8,$
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_a2 = 	WIDGET_LABEL(base_t3d, VALUE=' a2 (kg/l) ')


;ok = SetWidget_maxsizes([label_frame_period, label_nframes, label_frame_ini, label_dose, label_tr, label_r1], /XSIZE)
;ok = SetWidget_maxsizes([label_m1, label_m2, label_a1, label_a2,label_ktrans_ini,label_kep_ini], /XSIZE)

st_out ={$
	wbase_options    : base_y,$
	wtext_frame_period : text_frame_period,$
	wtext_frame_injection  : text_frame_injection,$
	wtext_nframes     : text_nframes,  $
	wtext_nframes_auc : text_nframes_auc,  $
	wtext_dose        : text_dose,$
	wtext_t10_tissue  : text_t10_tissue,$
	wtext_t10_blood   : text_t10_blood,$
	wtext_haematocrit : text_haematocrit,$
	wtext_flip_angle  : text_flip_angle,$
	wtext_r1          : text_r1,$
	wtext_tr          : text_tr,$
	wtext_m1          : text_m1,$
	wtext_m2          : text_m2,$
	wtext_a1          : text_a1,$
	wtext_a2          : text_a2,$
	wtext_ktrans_rr   : text_ktrans_rr,$
	wtext_kep_rr      : text_kep_rr $
}

ok= fwdcemri_Parameters_write(ST_WIDGETS=st_out, ST_PARAM=st_param)
ok= fwdcemri_ParametricAIF_write(ST_WIDGETS=st_out, ST_AIF=st_aif)

RETURN, st_out

END

;********************************************************************************************************
;********************************************************************************************************

FUNCTION Interface_DCEMRI_about, GROUP_LEADER=group_leader
	;-------------------------------------------------------------------
	file_image_logo  = '.\Icons\Logobit_small.png'
	IF FILE_TEST(file_image_logo) THEN BEGIN
		image_logo = READ_PNG(file_image_logo)
	ENDIF
	;-------------------------------------------------------------------
	strarr_labelinfo = STRARR(1,7)
	strarr_labelinfo[0] = 'DCE@urLAB v1.0'
	strarr_labelinfo[1] = 'T1-weighted DCE-MRI Pharmacokinetic analysis, version 1.1a, 30/01/2013'
	strarr_labelinfo[3:5] = Get_BITUPM_strarray()
	;strarr_labelinfo[3:5] = Get_License_Array()
	arr_font_type='arial'
    arr_font_effect=['bold','light','light','light','light','light','light']
	arr_font_size=[16,14,10,16,16,16,14]
	;arr_idl_info = [[Get_OS_info()], [Get_IDL_License_Info()]]
	arr_idl_info = [[Get_OS_info()], ['']]
	n_info = N_ELEMENTS(arr_idl_info)
	strarr_labelinfo= [[strarr_labelinfo],[arr_idl_info]]
	arr_font_effect = [arr_font_effect, REPLICATE('bold', n_info)]
	arr_font_size   = [arr_font_size, REPLICATE(12, n_info)]
	ok = Interface_LabelInfo(GROUP_LEADER=group_leader, TITLE='about', LABEL=strarr_labelinfo, $
		FONT_TYPE=arr_font_type, FONT_SIZE=arr_font_size, MODAL=1,$
		FONT_EFFECT=arr_font_effect, IMAGE_LOGO=image_logo)

	RETURN, 1
END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_DCEMRI_killnotify

	ok = DIALOG_MESSAGE('Quit DCEurLAB?', /QUESTION)
	IF STRUPCASE(ok) EQ 'NO' THEN RETURN
END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_DCEMRI_destroy, ev


HEAP_GC

WIDGET_CONTROL, ev, GET_UVALUE=pst

IF (*Pst).opt_nodestroyatend EQ 0 THEN BEGIN
	PTR_FREE, (*Pst).data
	PTR_FREE, (*Pst).mask
	PTR_FREE, (*Pst).results.RCE_im
	PTR_FREE, (*Pst).t1map
	PTR_FREE, (*Pst).results.dataptr
	PTR_FREE, (*Pst).results.arr_Cp
	PTR_FREE, (*Pst).results.arr_Signal_AIF
	PTR_FREE, (*Pst).results.arr_signal_RR
	PTR_FREE, (*Pst).results.arr_time
	PTR_FREE, (*Pst).results.data_CT
	PTR_FREE, (*Pst).results.arr_Ct_RR
	PTR_FREE, (*Pst).results.data_Signal_lr
	PTR_FREE, (*Pst).results.data_Signal_ROI
	PTR_FREE, (*Pst).roi.mask_lr
	PTR_FREE, (*Pst).roi.idx_points
	PTR_FREE, (*Pst).roi.idx_points_lr
	PTR_FREE, (*Pst).roi.im_idxs
	PTR_FREE, (*Pst).inf.slope
	PTR_FREE, (*Pst).inf.info
ENDIF

IF OBJ_VALID((*Pst).vi1.ocontainer)  THEN BEGIN
    OBJ_DESTROY, (*Pst).vi1.ocontainer
ENDIF
IF OBJ_VALID((*Pst).vi2.ocontainer)  THEN BEGIN
	OBJ_DESTROY, (*Pst).vi2.ocontainer
ENDIF

HEAP_GC
WIDGET_CONTROL,  (*Pst).wd.wbase_main,  /DESTROY
HEAP_GC
LOADCT, 0
;HELP, /HEAP

END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_DCEMRI_event, ev

opt_catcherror = 0
;*****************************************
IF LMGR(/RUNTIME) OR LMGR(/VM) OR (opt_catcherror EQ 1) THEN BEGIN
	CATCH, theError
	IF theError NE 0 THEN BEGIN
	  CATCH, /Cancel
	  void = DIALOG_MESSAGE(!ERROR_STATE.MSG, /ERROR)
	  RETURN
	ENDIF
ENDIF
;*****************************************

IF (TAG_NAMES(ev, /STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST') THEN BEGIN
	WIDGET_CONTROL, sEvent.top, /DESTROY
	RETURN
ENDIF

WIDGET_CONTROL, ev.top, GET_UVALUE=Pst
WIDGET_CONTROL, ev.id,  GET_UVALUE=uval
HEAP_GC

IF N_ELEMENTS(uval) EQ 0 THEN RETURN

;---------------------------;
; only to breakpoints, comment..
;IF ev.type NE 2 THEN BEGIN
;	PRINT, 'Breakpoint'
;ENDIF
;---------------------------
CASE uval[0] OF
'CLOSE_WINDOWS' : BEGIN
	WDEL_ALL
	RETURN
END
'BOPEN_MRI_DT1'  : BEGIN ; Abre ficheros

    WIDGET_CONTROL, /HOURGLASS
    uname = WIDGET_INFO(ev.id, /UNAME)

	size_old = (*Pst).inf.size_im ; backup to avoid reset the ROI if possible

    CASE uname OF
    'OPEN_MRI_RAW': BEGIN

        str_file = DIALOG_PICKFILE(DISPLAY_NAME='Open MRI-T1 dynamic data', DIALOG_PARENT=ev.top, $
         	/MUST_EXIST, PATH=(*Pst).path_data, TITLE = 'Open dynamic data')
    	IF FILE_TEST(str_file) EQ 0 THEN RETURN

		IF STRUPCASE(get_name_extension(str_file)) EQ 'RAW' THEN opt_inv_y=0 ELSE opt_inv_y=1

		st_data = {$
			xsize:(*Pst).inf.size_im[0],ysize:(*Pst).inf.size_im[1],$
			nslices:(*Pst).inf.nslices, nframes:(*Pst).par.nframes, offset:0, $
			frame_period:(*Pst).par.frame_period*60, frame_injection:(*Pst).par.nframe_injection,$
			typedata:12, invarr : [0,opt_inv_y,0,0l], littleendian:0}
		;---------------------------------------------------------------------------------
		std = Interface_ReadRaw(GROUP_LEADER=ev.top,FILE=str_file, ST_DATA=st_data)
		;---------------------------------------------------------------------------------
		IF SIZE(std, /TNAME) NE 'STRUCT' THEN RETURN

		WIDGET_CONTROL, /HOURGLASS
		;---------------------------------------------------------------------------------
		data = Read_raw_v2(FILE=str_file, OFFSET=std.offset, $
			DIMENSIONS=[std.xsize, std.ysize, std.nslices, std.nframes],$
			SWAP=std.littleendian, TYPE=std.typedata, REVERSE_ARR=std.invarr,$
			ERROR_ID = error_id, ERROR_STR = error_str,$
			MATCH_SIZE=0)
		IF SIZE(data, /N_DIMENSIONS) LE 2 THEN BEGIN
	    	ok = DIALOG_MESSAGE('Error reading data', /INFO) & 	RETURN
		ENDIF
		;---------------------------------------------------------------------------------
		(*Pst).inf.file         = str_file
		(*Pst).inf.size_im[0]   = std.xsize
		(*Pst).inf.size_im[1]   = std.ysize
		(*Pst).inf.nslices 		= std.nslices
		(*Pst).inf.nframes 		= std.nframes
		*(*Pst).inf.slope  		= 1
		(*Pst).par.frame_period = std.frame_period/60d ; en minutos
		(*Pst).par.nframe_injection = std.frame_injection
		*(*Pst).inf.info = ['Raw file', (*Pst).inf.file]

		data = FLOAT(data)*(*(*Pst).inf.slope) ; before load image...importante
	END
	'OPEN_MRI_BRUKER' : BEGIN
		str_path = DIALOG_PICKFILE(DISPLAY_NAME='Open MRI-T1 (Bruker Directory)', DIALOG_PARENT=ev.top, $
         	/MUST_EXIST, /DIRECTORY, PATH=(*Pst).path_data, TITLE = 'Open Bruker directory of dynamic data')
    	IF FILE_TEST(str_path, /DIRECTORY) EQ 0 THEN RETURN

		;-------------------------------------------------------------------------------
		data = Read_BrukerBiospin_image(str_path, ST_DATA=st_data, ERROR_ID=error_id, $
			ERROR_STR=error_str, ONLY_HEADER=opt_only_header, NOSLOPE=0)
		;-------------------------------------------------------------------------------
		IF error_id NE 0 THEN BEGIN
			ok = DIALOG_MESSAGE([['Error reading data'], [error_str]], /INFO) &	RETURN
		ENDIF
		IF st_data.nfrec NE st_data.nphases THEN BEGIN
			dim_max = MAX([st_data.nfrec, st_data.nphases], MIN=dim_min)
			str_dialog = [['Not square slices [' + STRTRIM(st_data.nfrec,2)+','+ STRTRIM(st_data.nphases,2) + ']'],$
				[' Convert to ' + STRTRIM(dim_max,2) +','+ STRTRIM(dim_max,2) + '] ?']]
			ok = DIALOG_MESSAGE(str_dialog, /QUESTION)
			IF STRUPCASE(ok) EQ 'YES' THEN BEGIN
				IF dim_max MOD dim_min EQ 0 THEN BEGIN
					data = REBIN(data,   dim_max,dim_max, st_data.nslices, st_data.nframes, /SAMPLE)
				ENDIF ELSE BEGIN
					data = CONGRID(data, dim_max,dim_max, st_data.nslices, st_data.nframes, /INTERP)
				ENDELSE
				st_data.nfrec   = dim_max
				st_data.nphases = dim_max
			ENDIF
		ENDIF
		;---------------------------------
		(*Pst).inf.size_im[0]   = st_data.nfrec
		(*Pst).inf.size_im[1]   = st_data.nphases
		(*Pst).inf.nslices    = st_data.nslices
		(*Pst).inf.nframes    = st_data.nframes
		*(*Pst).inf.slope     = st_data.visu_dataslope
		(*Pst).inf.file       = st_data.file
		(*Pst).par.frame_period = st_data.frametime/60 ; en minutos
		*(*Pst).inf.info      = st_data.info
		;----------------------------------------------------------------------------------------------
		array_questions = ['Injection frame (from 1 to n)']
		initial_result  = [STRTRIM(STRING((*Pst).par.nframe_injection+1),2)]
		result = Interface_datain_generic(GROUP_LEADER=ev.top, XSIZE=xsize, YSIZE=ysize, POSITION=position, $
	   		TITLE='Parameters', ARRAY_QUESTIONS=array_questions, INITIAL_RESULT=initial_result, ROW=0)
		;----------------------------------------------------------------------------------------------
		IF result[0] EQ '' THEN RETURN
		(*Pst).par.nframe_injection = (FIX(result[0])-1) > 0
	END
	'OPEN_MRI_DICOM' : BEGIN

		;PRINT, !VERSION.RELEASE
		;PRINT, !VERSION.OS
		IF !VERSION.MEMORY_BITS NE 32 THEN BEGIN
		ok = DIALOG_MESSAGE('DICOM images only can be read in the 32 bits version', /INFO)
			RETURN
		ENDIF

		str_path = DIALOG_PICKFILE(DISPLAY_NAME='Open DICOM dynamic data', DIALOG_PARENT=ev.top, $
         	/MUST_EXIST, PATH=(*Pst).path_data, TITLE = 'Open DICOM dynamic data', /DIRECTORY)
    	IF FILE_TEST(str_path, /DIRECTORY) EQ 0 THEN RETURN

    	WIDGET_CONTROL, /HOURGLASS
    	  	;----------------------------------------------------------------------------------
		data = Read_DynamicDicom(PATH=str_path, ERROR_STR=error_Str, ERROR_ID=error_id)
		;----------------------------------------------------------------------------------
		IF error_id NE 0 THEN BEGIN
			ok = DIALOG_MESSAGE([['Error reading data'], [error_str]], /INFO) &	RETURN
		ENDIF
		IF SIZE(data, /N_DIMENSIONS) NE 4 THEN BEGIN
			ok = DIALOG_MESSAGE('Error in data dimensions', /INFO) & RETURN
		ENDIF

		;---------------------------------------------------------------------------------
		array_questions = ['Injection frame (from 1 to n)', 'Frame period (seconds)']
		initial_result  = [STRTRIM(STRING((*Pst).par.nframe_injection+1),2), STRTRIM(STRING((*Pst).par.frame_period*60),2)]
		result = Interface_datain_generic(GROUP_LEADER=ev.top, XSIZE=xsize, YSIZE=ysize, POSITION=position, $
	   		TITLE='Parameters', ARRAY_QUESTIONS=array_questions, INITIAL_RESULT=initial_result, ROW=0)
		;---------------------------------------------------------------------------------
		IF result[0] EQ '' THEN RETURN
		(*Pst).par.nframe_injection = (FIX(result[0])-1) > 0
		(*Pst).par.frame_period =     DOUBLE(result[1])/60.0 ; en minutos

		;---------------------------------------------------------------------------------
		;str_1 = get_lastpath(get_path(STRMID(str_path, 0, STRLEN(str_path)-1)))
		;str_1 = STRMID(str_1, 0, STRLEN(str_1)-1)
		;HELP, data
		;ok = write_raw(FILE='F:\Experiments_DCE-MRI\QIBA (DCE-MRI test data sets)\QIBA_v07_Tofts\QIBA_v7_Tofts_Siemens\' + str_1 + '_(50x80x721)(uint).raw', DATA=data)
		;PRINT, OK & RETURN
		;---------------------------------------------------------------------------------
		(*Pst).inf.file         = str_path
		(*Pst).inf.size_im[0]   = (SIZE(data, /DIMENSIONS))[0]
		(*Pst).inf.size_im[1]   = (SIZE(data, /DIMENSIONS))[1]
		(*Pst).inf.nslices 		= (SIZE(data, /DIMENSIONS))[2]
		(*Pst).inf.nframes 		= (SIZE(data, /DIMENSIONS))[3]
		*(*Pst).inf.slope  		= 1
		*(*Pst).inf.info = ['DICOM files', (*Pst).inf.file]

		data = FLOAT(data)*(*(*Pst).inf.slope) ; before load image...importante

	END
	ELSE : RETURN
	ENDCASE

	ok= fwdcemri_Parameters_write(ST_WIDGETS=(*Pst).wdp, ST_PARAM=(*Pst).par)
	;-----------------------------------------------------------------------
	ok = fwdcemri_Load_image(Pst, data, FILENAME=(*Pst).inf.file, ERROR_STR=error_str)
    IF ok[0] NE 1 THEN BEGIN
        ok=DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main, error_str, /ERROR)
        RETURN
    ENDIF

    ok = fwdcemri_activate_bases(Pst, /OPEN_DATA)
    (*Pst).path_data = get_path((*Pst).inf.file)

	;(*Pst).vi1.oimageROI -> setProperty, HIDE=1

    WIDGET_CONTROL, (*Pst).wd.wbase_main, UPDATE=0
    ok = fwdcemri_Reshape_Windows(Pst, /TWO_WINDOWS)

    ok = fwdcemri_Draw_Image(Pst, CONTROL_ZOOM=(*Pst).initial_zoom, NODRAW=1)
    ok = fwdcemri_Draw_RCEImage(Pst, CONTROL_ZOOM=(*Pst).initial_zoom, /CALCULATE, NODRAW=1)

	ok = fwdcemri_DrawInSameSize(Pst)
    ok = fwdcemri_Redraw(Pst)
    ;ok = fwdcemri_checkDrawRoi(Pst)
    WIDGET_CONTROL, (*Pst).wd.wbase_main, UPDATE=1

    (*Pst).vi1.Opolyline_ROI -> getProperty, DATA=datapolyline
    IF N_ELEMENTS(datapolyline) GT 2 THEN BEGIN
	    IF TOTAL(size_old EQ (*Pst).inf.size_im, /INTEGER) NE 2 THEN BEGIN
	    	OK = DIALOG_MESSAGE('Different size to the previous data: ROI reset', /INFO)
	    	(*Pst).flag_noresetroi = 0
	    	WIDGET_CONTROL, (*Pst).wd.wbttn_newROI, SEND_EVENT={ID:0L, TOP:0L, HANDLER:0L}
		ENDIF
	ENDIF
END
;-------------------------------------------------------------------------------------------------
'IMAGE_INFO'   : BEGIN

	IF N_ELEMENTS(*(*Pst).data) LE 2 THEN BEGIN
		ok = DIALOG_MESSAGE('No image loaded', /INFO) & RETURN
	ENDIF

	strarr_info = fwdcemri_ImageInfo(Pst, TITLE=title)
	;ok = DIALOG_MESSAGE(REFORM(strarr_info), /INFO, /CENTER, TITLE=title)
	XDISPLAYFILE, TEXT=strarr_info, TITLE='Information', DONE_BUTTON='ok', GROUP=ev.top, HEIGHT=N_ELEMENTS(strarr_info)+2, RETURN_ID=widget_id_info
	ok= GetWidget_PositionInScreen(widget_id_info, CENTER=1)


END
;-------------------------------------------------------------------------------------------------
'BOPEN_T1MAP'  : BEGIN ; Abre ficheros de mapa T1

	WIDGET_CONTROL, /HOURGLASS
	str_file = DIALOG_PICKFILE(DISPLAY_NAME='Open T1 map', DIALOG_PARENT=ev.top, $
		/MUST_EXIST, PATH=(*Pst).path_t1map, TITLE = 'Open T1 map')

	str_extension = get_name_extension(str_file)
	CASE STRUPCASE(str_extension) OF
	'RAW' : BEGIN
		;'F:\Experiments_DCE-MRI\dMRI_rui090928.Vu6\Results\';
		t1map = read_raw(FILE=str_file, DIMENSIONS=[(*Pst).inf.size_im[0], (*Pst).inf.size_im[1], (*Pst).inf.nslices], TYPE='FLOAT')

	END
	'' : BEGIN ; formato de bruker procesado, leyendo la cabecera para el slope

		file_header = get_path(str_file) + 'visu_pars'
		IF FILE_TEST(file_header) NE 1 THEN BEGIN
			ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main,'Header "visu_pars" not found in folder', /INFO)
			RETURN
		ENDIF
		st_visuheader = Read_BrukerBiospin_VisuHeader(file_header, ERROR_ID=error_id)
		IF error_id NE 0 THEN BEGIN
			ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main,'Error reading header, not compatible format', /INFO)
			RETURN
		ENDIF
		dataslope = st_visuheader.dataslope

		nimages = 5
		IF N_ELEMENTS(dataslope) NE nimages*(*Pst).inf.nslices THEN BEGIN
			ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main,'t1 map not compatible with image size', /INFO)
			RETURN
		ENDIF

		dataslope_sel = dataslope[INDGEN((*Pst).inf.nslices)*nimages+2] ; EL DOS DICE QUE LA IMAGEN ES LA Nº 3
		t1map = Read_bruker_T1map(str_file, ROTATION=1, DIMENSIONS=[(*Pst).inf.size_im[0], (*Pst).inf.size_im[1], (*Pst).inf.nslices], $
			DATA_SLOPE=dataslope_sel, DATA_OFFSET=data_offset, ERROR_STR=error_str)

	END
	ELSE : BEGIN
		ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main,'format not supported', /INFO)
		RETURN
	END
	ENDCASE

	IF (SIZE(t1map, /N_DIMENSIONS) NE 3) THEN BEGIN
		ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main,'Error reading T1 map file, dimensions does not agree with data', /INFO)
		RETURN
	ENDIF
	PRINT, MAX(t1map)
	PRINT, MIN(t1map)

	(*Pst).path_t1map = get_path(str_file)
	*(*Pst).t1map = t1map/1000.0 ; pasamos el tiempo de relajación (t1) a segundos !!

	ok = fwdcemri_activate_bases(Pst, /OPEN_T1MAP)

	;ok = fwdcemri_checkDrawRoi(Pst)

END
;-------------------------------------------------------------------------------------------------
'OPEN_CURVE_TXT'  : BEGIN


	uname = WIDGET_INFO(ev.id, /UNAME)

	str_file = DIALOG_PICKFILE(DISPLAY_NAME='Open curve file', DIALOG_PARENT=ev.top, $
	   	/MUST_EXIST, PATH=(*Pst).path_aif, TITLE = 'Open curve file')
	IF FILE_TEST(str_file) EQ 0 THEN RETURN
	(*Pst).path_aif = get_path(str_file)
	arr_Signal      = Read_AIFfunction(str_file, ARR_TIME=arr_time, ARR_SIGNAL=arr_data, ERROR_ID=error_id, ERROR_STR=error_str)

	IF error_id NE 0 THEN BEGIN & ok = DIALOG_MESSAGE(error_str, /INFO) & RETURN &	ENDIF

	; PENDIENTE: preguntar si se quiere ajustar el numero de frames y el periodo a la señal.
	IF N_ELEMENTS(arr_Data) LT (*Pst).inf.nframes THEN BEGIN
		ok = DIALOG_MESSAGE('Dynamic Curve contains less values than image frames. Replicate last value?', /QUESTION)
		IF STRUPCASE(ok) EQ 'YES' THEN BEGIN
			last_value = arr_data[N_ELEMENTS(arr_Data)-1]
			arr_data = [arr_data, REPLICATE(last_value,  (*Pst).inf.nframes-N_ELEMENTS(arr_Data))]
		ENDIF ELSE BEGIN
			ok = DIALOG_MESSAGE('Curve not loaded', /INFO) &	RETURN
		ENDELSE
	ENDIF
	IF N_ELEMENTS(arr_Data) GT (*Pst).inf.nframes THEN BEGIN
		ok = DIALOG_MESSAGE('Array of data contains more values than time samples (Loading first values)', /INFO)
	ENDIF

	CASE uname OF
    'AIF_CURVE' : BEGIN
		*(*Pst).results.arr_Signal_AIF = arr_data;                  [0:(*Pst).inf.nframes-1]
		*(*Pst).results.arr_time_AIF   = arr_time/60d; ; en minutos [0:(*Pst).inf.nframes-1]
		(*Pst).flag_AIF_from_file      = 1
		(*Pst).flag_AIF_from_ROI       = 0
		method_plot =3
	END
    'RR_CURVE'  : BEGIN
		*(*Pst).results.arr_Signal_RR  = arr_data;                  [0:(*Pst).inf.nframes-1]
		*(*Pst).results.arr_time_RR    = arr_time/60d; ; en minutos [0:(*Pst).inf.nframes-1]
		(*Pst).flag_RR_from_file      = 1
		(*Pst).flag_RR_from_ROI       = 0
		method_plot = 3
	END
	ENDCASE
	win_number = get_windowNumber(OPTNEXT=1, XPOS=xpos, YPOS=ypos)
	ok = Plot_profile_dynMRI(arr_data, ARR_TIME=arr_time/60d, WIN_NUMBER=win_number, METHOD=method_plot)

END
;-------------------------------------------------------------------------------------------------

'OPEN_CURVE_PARAMETRIC'  : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)
    CASE uname OF
	'AIF_PARAMETRIC' : BEGIN
		str_file = DIALOG_PICKFILE(DISPLAY_NAME='Open parametric AIF config file', DIALOG_PARENT=ev.top, $
		   	/MUST_EXIST, PATH=(*Pst).path_aif, TITLE = 'Open parametric AIF config file')
	    IF FILE_TEST(str_file) EQ 0 THEN RETURN

		(*Pst).path_aif = get_path(str_file)
		st_cpmodel = Read_CpModelConfig(str_file, MODEL=model, ERROR_ID=error_id, ERROR_STR=error_str)
		IF error_id NE 0 THEN BEGIN
			ok = DIALOG_MESSAGE(error_str, /INFO) & RETURN
		ENDIF

		IF st_cpmodel.model NE 'biexponential' THEN BEGIN
			ok = DIALOG_MESSAGE('Only biexponential model allowed in this version', /INFO) & RETURN
		ENDIF

		*(*Pst).parametric_AIF = TEMPORARY(st_cpmodel) ; Store the structure in the pointer reserved to the parametric AIF

		ok= fwdcemri_ParametricAIF_write(ST_WIDGETS=(*Pst).wdp, ST_AIF=*(*Pst).parametric_AIF)

		win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)

		arr_Cp = ParametricCp_to_CurveCp(ST_AIF=*(*Pst).parametric_AIF, ARR_TIME=*(*Pst).results.arr_time, DOSE=(*Pst).par.dose, INJECTION_FRAME=(*Pst).par.nframe_injection)
		ok = Plot_profile_dynMRI(arr_Cp, ARR_TIME=*(*Pst).results.arr_time, WIN_NUMBER=win_number, METHOD=3)
	END
	ENDCASE
END
;-------------------------------------------------------------------------------------------------
'CLOSE_AIF_TXT' : BEGIN
		(*Pst).flag_AIF_from_ROI       = 0
		(*Pst).flag_AIF_from_file      = 0
		*(*Pst).results.arr_Signal_AIF = -1
		*(*Pst).results.arr_time_AIF   = -1  ; destroy AIF curves
END
;-------------------------------------------------------------------------------------------------
'SELECT_ROI' : BEGIN ; Open the AIF interface

	;inc_time = (*(*Pst).results.arr_time)[1]-(*(*Pst).results.arr_time)[0]

	IF (*Pst).flag_drawROI NE 1 THEN BEGIN
		ok = DIALOG_MESSAGE('No ROI selected', /INFO)
		RETURN
	ENDIF
	arr_dyn  = fwdcemri_DynDataFromIdx(Pst)
	IF arr_dyn[0] EQ -1 THEN BEGIN
		ok = DIALOG_MESSAGE('No valid or selected ROI', /INFO)
		RETURN
	ENDIF
	arr_dyn_mean = TOTAL(arr_dyn,1, /DOUBLE)/(SIZE(arr_dyn,/DIMENSIONS))[0]

	uname= WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'AIF_FROM_ROI' : BEGIN

		*(*Pst).results.arr_Signal_AIF = TEMPORARY(arr_dyn_mean)
		*(*Pst).results.arr_time_AIF   = *(*Pst).results.arr_time
		(*Pst).flag_AIF_from_file      = 0
		(*Pst).flag_AIF_from_ROI       = 1
		ok = DIALOG_MESSAGE('AIF curve loaded from actual ROI', /INFO)

	END
	'RR_FROM_ROI' : BEGIN
		*(*Pst).results.arr_signal_RR  = TEMPORARY(arr_dyn_mean)
		*(*Pst).results.arr_time_RR   = *(*Pst).results.arr_time
		(*Pst).flag_RR_from_file      = 0
		(*Pst).flag_RR_from_ROI       = 1
		ok = DIALOG_MESSAGE('RR curve loaded from actual ROI', /INFO)
	END
	ENDCASE

END
;-------------------------------------------------------------------------------------------------
'BFREE_DATA' : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)

    CASE uname[0] OF
    'FREE_MRI_T1': BEGIN
    	HELP, /MEMORY
    	(*Pst).vi1.oimage->SetProperty, DATA=BYTARR(10,10)
		(*Pst).vi2.oimage->SetProperty, DATA=BYTARR(10,10)
		(*Pst).flag_drawroi = 0
       	ok = fwdcemri_Redraw(Pst)
       	ok = fwdcemri_Hide_Images(Pst)
       	undefine, *(*Pst).data
       	ok = fwdcemri_ResetResults(Pst)
       	ok = fwdcemri_activate_bases(Pst, /FREE_DATA)
    	HELP, /MEMORY
    END
	'FREE_T1MAP': BEGIN
		undefine, *(*Pst).t1map
		ok = fwdcemri_activate_bases(Pst, /FREE_T1MAP)
		;ok = fwdcemri_checkDrawRoi(Pst)
	END
	ENDCASE
END
;-------------------------------------------------------------------------------------------------------
'CHOOSE_MODEL' : BEGIN

	st_par = (*Pst).par
	ok= fwdcemri_Parameters_read (ST_WIDGETS=(*Pst).wdp, ST_PARAM=st_par, ERROR_ID=error_id, ERROR_STR=error_str)
	ok= fwdcemri_Parameters_write(ST_WIDGETS=(*Pst).wdp, ST_PARAM=st_par)
	IF error_id NE 0 THEN BEGIN
		ok = DIALOG_MESSAGE(error_str, /INFO) &	RETURN ; no es una estructura
	ENDIF
	(*Pst).par = st_par
	*(*Pst).results.arr_time = INDGEN((*Pst).par.nframes)*(*Pst).par.frame_period

END
;-------------------------------------------------------------------------------------------------------
'CHOOSE_PLASMACURVE' : BEGIN

	;strarr_plasmacurve = ['Not selected', 'Analytic biexponential', 'From File', 'From Arterial ROI']
	opt_dlist = WIDGET_INFO((*Pst).wd.wdropls_AIFcurve, /DROPLIST_SELECT)
	error_select = 0l

	CASE opt_dlist OF
	0 : (*Pst).opt_AIFcurve = opt_dlist
	1 : BEGIN
	 	IF SIZE(*(*Pst).parametric_aif, /TNAME) EQ 'STRUCT' THEN BEGIN
	 		IF (*(*Pst).parametric_aif).model EQ 'biexponential' THEN BEGIN
	 			ok = DIALOG_MESSAGE('Using AIF biexponential model:, parameters in AIF tab', /INFO)
	 			(*Pst).opt_AIFcurve = opt_dlist
	 		ENDIF ELSE BEGIN
	 			ok = DIALOG_MESSAGE('Error, non-exponential model is not supported', /INFO)
	 			error_select = 1
	 		ENDELSE
	 	ENDIF ELSE BEGIN
	 		ok = DIALOG_MESSAGE('No biexponential AIF model has been loaded', /INFO)
	 		error_select     = 2
		ENDELSE
	END
	2 : BEGIN
		IF (N_ELEMENTS(*(*Pst).results.arr_signal_AIF) LE 1) OR ((*Pst).flag_aif_from_file EQ 0) THEN BEGIN
			ok = DIALOG_MESSAGE('No AIF curve has been loaded from file', /INFO)
		   	error_select = 3
		ENDIF ELSE BEGIN
			(*Pst).opt_AIFcurve = opt_dlist
		ENDELSE
	END
	3 : BEGIN
		IF (N_ELEMENTS(*(*Pst).results.arr_signal_AIF) LE 1) OR ((*Pst).flag_aif_from_ROI EQ 0) THEN BEGIN
			ok = DIALOG_MESSAGE('No AIF curve has been loaded from ROI', /INFO)
		   error_select = 4
		ENDIF ELSE BEGIN
			(*Pst).opt_AIFcurve = opt_dlist
		ENDELSE
	END
	ENDCASE

	IF error_select NE 0 THEN BEGIN
		WIDGET_CONTROL, (*Pst).wd.wdropls_AIFcurve, SET_DROPLIST_SELECT=0
	 	(*Pst).opt_AIFcurve = 0
	ENDIF
	WIDGET_CONTROL, (*Pst).wd.wbttn_viewAIF, SENSITIVE=error_select EQ 0


END
;-------------------------------------------------------------------------------------------------------
'CHOOSE_RRCURVE' : BEGIN

	;strarr_rrcurve = ['Not selected',     'Loaded from disk', 'From RR ROI']
	opt_dlist = WIDGET_INFO((*Pst).wd.wdropls_RRcurve, /DROPLIST_SELECT)

	error_select = 0l

	CASE opt_dlist OF
	0 : (*Pst).opt_RRcurve = opt_dlist
	1 : BEGIN
		IF (N_ELEMENTS(*(*Pst).results.arr_signal_RR) LE 1) OR ((*Pst).flag_RR_from_file EQ 0) THEN BEGIN
			ok = DIALOG_MESSAGE('No RR curve has been loaded from file', /INFO)
		   	error_select = 3
		ENDIF ELSE BEGIN
			(*Pst).opt_RRcurve = opt_dlist
		ENDELSE
	END
	2 : BEGIN
		IF (N_ELEMENTS(*(*Pst).results.arr_signal_RR) LE 1) OR ((*Pst).flag_RR_from_ROI EQ 0) THEN BEGIN
			ok = DIALOG_MESSAGE('No RR curve has been loaded from ROI', /INFO)
		   error_select = 4
		ENDIF ELSE BEGIN
			(*Pst).opt_RRcurve = opt_dlist
		ENDELSE
	END
	ENDCASE

	IF error_select NE 0 THEN BEGIN
		WIDGET_CONTROL, (*Pst).wd.wdropls_RRcurve, SET_DROPLIST_SELECT=0
	 	(*Pst).opt_RRcurve = 0
	ENDIF
	WIDGET_CONTROL, (*Pst).wd.wbttn_viewRR, SENSITIVE=error_select EQ 0

END
;-------------------------------------------------------------------------------------------------------
'EXPORT_ASCII' : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)

	reverse_lines    = 1  ; por defecto
	opt_index_points = 1L ; opción de poner los indices de los puntos
	opt_name         = 0
	str_info = ''

	CASE uname OF
	'RCE_VALUES' : BEGIN  ; cinética temporal
    	win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
		im_temp = fwdcemri_ROIkinetics(Pst, WIN_NUMBER=win_number,XPOS=xpos, YPOS=ypos, TYPE_PLOT=3);
		str_file_ini = 'ROI_kinetics_RCE'
    	reverse_lines    =  0 ;
    	opt_index_points = 0l
		opt_name      = 3
		str_label = '##  Time (s) ' + STRING(9B) + 'Average (AU) ' + STRING(9B) + 'SD'
		undefine, identifier
	END
	'ABS_VALUES': BEGIN
		win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
		im_temp = fwdcemri_ROIkinetics(Pst, WIN_NUMBER=win_number,XPOS=xpos, YPOS=ypos, TYPE_PLOT=3, /ABSOLUTE_VALUE);
		str_file_ini = 'ROI_kinetics_abs'
    	reverse_lines    =  0 ;
    	opt_index_points = 0l
		opt_name      = 3
		str_label = '##  Time (s) ' + STRING(9B) + 'Average (AU) ' + STRING(9B) + 'SD'
		undefine, identifier
	END
	ELSE:
    ENDCASE
    ;------------------------------------
	str_array = Translate_data2strings(im_temp, N_DECIMALS=5, REVERSE_LINES=reverse_lines, INTEGER=opt_integer)

	IF (SIZE(str_array,/TYPE) NE 7) THEN BEGIN
		IF str_array[0] EQ -1 THEN BEGIN
			ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main, 'Nothing to export', /INFO)
			RETURN
		ENDIF
	ENDIF
	IF N_ELEMENTS(str_label) NE 0 THEN str_array = [str_label, str_array]

	str_pointini = fwdcemri_suffixresultname(ST_ROI=(*Pst).roi, SIZEIMAGE=(*Pst).inf.size_im, '[ROIf]')

	str_file_ini+=str_pointini +'.txt

	str_file = DIALOG_PICKFILE(PATH=(*Pst).path_result, TITLE='Export ' + str_info  + ' (ASCII file)', $
		DEFAULT_EXTENSION = 'txt', WRITE=1, FILE=str_file_ini, DIALOG_PARENT=ev.top)
	IF str_file EQ '' THEN RETURN

	(*Pst).path_result = get_path(str_file)
	file_name  = get_path(str_file) + get_name_field(str_file) + '.txt'
	ok = WriteFile_ASCII(file_name, ARRAY=str_array, ERROR_ID=error_id)

	IF ok EQ 1 THEN BEGIN
		ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main,['Result saved in file ', get_name(file_name)], /INFO)
	ENDIF
END
;-------------------------------------------------------------------------------------------------
'INITIAL_PARAMETERS' : BEGIN
	PRINT, 'change initial parameters'

	st_ini_result = {ktrans_ini: (*Pst).par.ktrans_ini,$
		   kep_ini :  (*Pst).par.kep_ini, $
		   Ah_ini  :  (*Pst).par.Ah_ini,  $
		   keph_ini:  (*Pst).par.keph_ini,$
		   kelh_ini:  (*Pst).par.kelh_ini,$
		   sl_ini  :  (*Pst).par.sl_ini,  $
		   kepl_ini:  (*Pst).par.kepl_ini,$
		   ok : 0b $
	 } ; valores aquí en s-1, pero se representan en pantalla por min-1

	ok = Interface_InitialParameters(GROUP_LEADER=(*Pst).wd.wbase_main, POSITION=position, $
    	RESULT=result, INITIAL_RESULT=st_ini_result)

	IF ok EQ 1 THEN BEGIN
		(*Pst).par.ktrans_ini = result.ktrans_ini
		(*Pst).par.kep_ini    = result.kep_ini
		(*Pst).par.Ah_ini     = result.Ah_ini
		(*Pst).par.keph_ini   = result.keph_ini
		(*Pst).par.kelh_ini   = result.kelh_ini
		(*Pst).par.kepl_ini   = result.kepl_ini
		(*Pst).par.Sl_ini     = result.Sl_ini
	ENDIF
	PRINT, (*Pst).par.ktrans_ini,(*Pst).par.kep_ini, (*Pst).par.Ah_ini,(*Pst).par.keph_ini,(*Pst).par.kelh_ini,$
			(*Pst).par.Sl_ini, (*Pst).par.kepl_ini


END
;-------------------------------------------------------------------------------------------------
'DYNAMIC_RANGE' : BEGIN

	IF N_ELEMENTS(*(*Pst).data) LE 2  THEN RETURN

	;ok = DIALOG_MESSAGE('Unavailable option in this version', DIALOG_PARENT = (*Pst).wd.wbase_main, /INFO)
	;RETURN

	str_lastframe1   = STRTRIM(STRING((*Pst).inf.nframes-1),2)
	str_lastframe2   = STRTRIM(STRING((*Pst).inf.nframes),2)
	array_questions = ['Initial frame (1-'+str_lastframe1+')', 'Last frame (2-'+str_lastframe2+')']
	initial_result  = ['1', str_lastframe2]

	result = Interface_datain_generic(GROUP_LEADER=ev.top, XSIZE=xsize, YSIZE=ysize, POSITION=position, $
    	TITLE='Dynamic range', ARRAY_QUESTIONS=array_questions, INITIAL_RESULT=initial_result, ROW=0)
    ;----------------------------------------------------------------------------------------------
   IF MAX(STRLEN(result)) EQ 0 THEN RETURN

	ok = DIALOG_MESSAGE('Warning, the operation cannot be undone, continue?', /QUESTION)
	IF STRUPCASE(ok) EQ 'NO' THEN RETURN

	frame_ini = FIX(result[0]-1)
	frame_end = FIX(result[1]-1)

	IF (frame_ini LT 0) OR (frame_ini GT ((*Pst).inf.nframes-2)) OR $
		(frame_end LT 1) OR (frame_end GT ((*Pst).inf.nframes-1)) OR $
		(frame_ini GE frame_end-1) THEN BEGIN
		ok=DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main, 'Invalid range', /INFO)
		RETURN
	ENDIF

	IF (frame_ini NE 0) OR (frame_end NE ((*Pst).inf.nframes-1)) THEN BEGIN
		*(*Pst).data = (*(*Pst).data)[*,*,*,frame_ini:frame_end]
		ok = fwdcemri_Load_image(Pst, *(*Pst).data, FILENAME=str_file, ERROR_STR=error_str)
	    IF ok[0] NE 1 THEN BEGIN
	        ok=DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main, error_str, /ERROR)
	        RETURN
	    ENDIF
	    ok = fwdcemri_activate_bases(Pst, /OPEN_DATA)

	    ok = fwdcemri_Draw_RCEImage(Pst, /CALCULATE, NODRAW=1)
		ok = fwdcemri_Draw_Image(Pst, NODRAW=1)
		ok = fwdcemri_reDraw(Pst)
	ENDIF
END
;-------------------------------------------------------------------------------------------------
'VIEW_MRI' : BEGIN
	uname= WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'VIEW_DYN_T1': BEGIN
		WIDGET_CONTROL, (*Pst).wd.wslide_n,   SENSITIVE=1
	END
	'VIEW_MRI_T2'    : BEGIN
		WIDGET_CONTROL, (*Pst).wd.wslide_n,   SENSITIVE=0
	END
	ENDCASE

	ok = fwdcemri_Draw_Image(Pst, CONTROL_ZOOM=(*Pst).initial_zoom)

END
;-------------------------------------------------------------------------------------------------
; Escribe gifs animados...
'WRITE_GIF' : BEGIN
n_slices = 6l
n_frames = 41l
FOR i=0,n_slices-1 DO BEGIN
	dym_im =BYTSCL(REFORM((*(*Pst).data)[*,*,i,*]))
	filename_gif = 'Gif_n' + STRTRIM(i+1,2) + '.gif'

	FOR j=0, n_frames-1 DO BEGIN
		image = REBIN(dym_im[*,*,j],256,256)
		WRITE_GIF, filename_gif, image,DELAY_TIME=10, CLOSE=j EQ (n_frames-1),/MULTIPLE, REPEAT_COUNT=0
	ENDFOR
ENDFOR
END

;-------------------------------------------------------------------------------------------------
'BVIEW_T1MAP' : BEGIN

	IF SIZE((*(*Pst).t1map), /N_DIMENSIONS) LT 2 THEN RETURN
	(*Pst).vi1.opalette->getproperty, BLUE_VALUES=b, GREEN_VALUES=g, RED_VALUES=r

	ok = view_mosaic (REFORM(*(*Pst).t1map), N_COL=3, DIV=0.5, WINDOW_NUMBER=0, XPOS=20, YPOS=20, $
		BACKGROUND='FFFFFF'x, BYTSCL_EQUAL=1, PALETTE=[[r],[g],[b]])
 END
;------------------------------------------------------------------------------------------------
'BVIEW_DYN_MRI' : BEGIN

	WIDGET_CONTROL,(*Pst).wd.wslide_z,  GET_VALUE=n_z
	nframes = (*Pst).inf.nframes

	IF SIZE((*(*Pst).data), /N_DIMENSIONS) LT 2 THEN RETURN

	uname= WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF

	'VIEW_DYN_MRI': BEGIN

		ok = DIALOG_MESSAGE('Unavailable option in this version', DIALOG_PARENT = (*Pst).wd.wbase_main, /INFO)
		RETURN

		n_col = FIX(CEIL(SQRT(nframes)))
		n_fil = (*Pst).inf.nframes/n_col
		DEVICE, Get_Screen_Size = screenSize
		div   = ((screensize[0]*1.0/n_fil)/(*Pst).inf.size_im[0])*1.3 > 0.5
		(*Pst).vi1.opalette->getproperty, BLUE_VALUES=b, GREEN_VALUES=g, RED_VALUES=r
		ok = view_mosaic (REFORM((*(*Pst).data)[*,*,n_z,*]), N_COL=n_col, DIV=div, WINDOW_NUMBER=0, $
			XPOS=20, YPOS=20, BACKGROUND='FFFFFF'x, BYTSCL_EQUAL=1, PALETTE=[[r],[g],[b]])
	END
	'VIEW_CINE' : BEGIN

		datacine = REFORM( (*(*Pst).data)[*,*,n_z,*])
		st_draw = WIDGET_INFO((*Pst).wd.wdraw1, /GEOMETRY)
		(*Pst).vi1.oimage->getproperty, DIMENSIONS=dim
		opalette =(*Pst).vi1.opalette

		Interface_Player, DATA=datacine, PALETTE=opalette, SIZE_VIEW=[st_draw.scr_xsize, st_draw.scr_ysize], SIZE_SCROLL=[st_draw.draw_xsize, st_draw.draw_ysize]


	END
	ENDCASE

	opt_write_sequence = 0
	IF opt_write_sequence THEN BEGIN
		max_abs = MAX((*(*Pst).data)[*,*,n_z,*])
		FOR i=0l, (*Pst).inf.nframes-1 DO BEGIN
			im = REFORM((*(*Pst).data)[*,*,n_z,i])
			imb = BYTSCL(im, MAX=max_abs)
			filename_png = 'E:\Download\' + 'Frame_' + Translate_number2string(NUMBER=i+1, LENGTH=3) + '.png'
			WRITE_PNG, Filename_png, imb
		ENDFOR
	ENDIF

END
;-------------------------------------------------------------------------------------------------
'S_CHANGE' : BEGIN   ; Elige la imagen 2D a mostrar (cambiando slice o frame)

    uname= WIDGET_INFO(ev.id, /UNAME)
    IF uname EQ 'SLICE_Z' THEN BEGIN
    	ok = fwdcemri_Draw_RCEImage(Pst)
	ENDIF
	ok = fwdcemri_Draw_Image(Pst)

	WIDGET_CONTROL, (*Pst).wd.wtext_ROI_size[0], SEND_EVENT={ID:0L, TOP:0L, HANDLER:0L}

	;---------------------------------------------------------------------------------------------
	; apply changes in another widget
	;pos_valid_widget = Locate_FirstValid_Widget((*Pst).wd.wAuxMainBases)
	;IF pos_Valid_widget NE -1 THEN BEGIN
	;	WIDGET_CONTROL, (*Pst).wd.wAuxMainBases[pos_valid_widget], GET_UVALUE=staux
	;	WIDGET_CONTROL, (*Pst).wd.wslide_n,    GET_VALUE=n_n
	;	WIDGET_CONTROL, staux.wd.wslide_n, SET_VALUE=n_n
	;	WIDGET_CONTROL, staux.wd.wslide_n, SEND_EVENT={ID:0L, TOP:staux.wd.wBase_main, HANDLER:staux.wd.wBase_main}
	;ENDIF
	;----------------------------------------------------------------------------------------------


END
;--------------------------------------------------------------------------------------------------
'VAL_CHANGE' : BEGIN

	uname= WIDGET_INFO(ev.id, /UNAME)

	IF uname EQ 'SLIDE_THRESHOLD' THEN opt_calculate = 1 ELSE opt_calculate = 0
	ok = fwdcemri_Draw_RCEImage(Pst, CALCULATE=opt_calculate)
END
;-------------------------------------------------------------------------------------------------------
'EXPORT_ROI' : BEGIN
	;ok = DIALOG_MESSAGE('You must perform ROI analysis and export from results interface', DIALOG_PARENT=(*Pst).wd.wbase_main, /INFO)


	roi_st = UpdateROIdata(IDX_POINTS=*(*Pst).roi.idx_points, ST_ROI=(*Pst).roi, IMAGESIZE=(*Pst).inf.size_im, SLICEZ=(*Pst).inf.current_slice, $
				ERROR_ID=error_id, ERROR_STR=error_str)
	IF error_id NE 0 THEN BEGIN
		ok = DIALOG_MESSAGE(error_str, /INFO) & RETURN & ENDIF
	(*Pst).roi = roi_st

	path_results = (*Pst).path_result
	ok = fwdcemri_ExportROI((*Pst).roi, IMAGE_SIZE = (*Pst).inf.size_im, PATH=path_results , OPOLYLINE= (*Pst).vi1.opolyline_ROI)
	(*Pst).path_result =path_results
END
;--------------------------------------------------------------------------------------------------
'IMPORT_ROI' : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'IMPORT_ROI_SAV' : BEGIN
		str_file   = DIALOG_PICKFILE(PATH=(*Pst).path_result, TITLE='import ROI (IDL-SAV format)', DIALOG_PARENT=ev.top)

		IF FILE_TEST(str_file) NE 1 THEN RETURN
		ok = fwdcemri_resetROI(Pst)
		RESTORE, str_file
		HELP, idx_points
		HELP, polylineROI
		;--------------------------------------------------------------
		; mantiene la compatibilidad con versiones anteriores
		IF N_ELEMENTS(ResolutionRoi) EQ 0 THEN ResolutionRoi = 1L
		IF N_ELEMENTS(typeROI) EQ 0 THEN typeROI = 2l
		;--------------------------------------------------------------
		HELP, resolutionROI
		HELP, typeROI
		*(*Pst).roi.idx_points = idx_points
		ok = fwdcemri_plotROI_set(Pst, DATA=polylineROI)
		(*Pst).roi.np_lr = ResolutionRoi
		(*Pst).roi.opt_typeROI = typeROI
		;--------------------------------------------------------------
		pos_resolution = (WHERE(STRTRIM(FIX(resolutionROI),2) EQ (*Pst).wd.wdropl_ROIresolution_array,ct))[0]
		IF ct NE 1 THEN BEGIN
			ok = DIALOG_MESSAGE('Saved resolution not supported: changed to 1', /INFO)
			pos_resolution = 0
		ENDIF
		WIDGET_CONTROL, (*Pst).wd.wdropl_ROIresolution, SET_DROPLIST_SELECT=pos_resolution
		;--------------------------------------------------------------
		WIDGET_CONTROL,(*Pst).wd.wdropl_typeROI, SET_DROPLIST_SELECT=(*Pst).roi.opt_typeROI

		CASE (*Pst).roi.opt_typeROI OF
		0 : BEGIN
			p_ini_x = *(*Pst).roi.idx_points[0] MOD (*Pst).inf.size_im[0]
			p_ini_y = *(*Pst).roi.idx_points[0] / (*Pst).inf.size_im[0]
			point_ini = fwdcemri_boxROI_set(Pst, POINT_INI=[p_ini_x,p_ini_y], DRAWTWO=1, $
		  		ID_TEXT=[(*Pst).wd.wtext_info1,(*Pst).wd.wtext_info2])
		END
		2 : BEGIN
			(*Pst).flag_freeROIclose      = 1
			(*Pst).flag_freeROIinit       = 0
			(*Pst).flag_freeROIcontinuous = 0
		END
		ELSE:
		ENDCASE
		(*Pst).path_result = get_path(str_file)
	END
	ELSE:RETURN
	ENDCASE

	ok = fwdcemri_activate_bases(Pst, /IMPORT_ROI)
	WIDGET_CONTROL, (*Pst).wd.wbttn_drawROI, SET_BUTTON=1

	;(*Pst).vi1.oimageROI -> setProperty, HIDE=1
	(*Pst).flag_noresetroi = 1

	ok = fwdcemri_redraw(Pst)

END
;-------------------------------------------------------------------------------------------------------
'ROI_KINETICS' : BEGIN

	uname= WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'TYPE_PLOT_1': BEGIN
		win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
		data_kinetics = fwdcemri_ROIkinetics(Pst, PLOT=1, WIN_NUMBER=win_number, XPOS=xpos, YPOS=ypos, TYPE_PLOT=1)
	END
	'TYPE_PLOT_3': BEGIN
		win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
		data_kinetics = fwdcemri_ROIkinetics(Pst, PLOT=1, WIN_NUMBER=win_number, XPOS=xpos, YPOS=ypos, TYPE_PLOT=3)
	END
	ENDCASE
END
;-------------------------------------------------------------------------------------------------------
'PLOT_AIF' : BEGIN

	CASE (*Pst).opt_AIFcurve OF
	0 : ; NOTHING
	1 : BEGIN ; analytic model for Cp
		IF SIZE(*(*Pst).parametric_AIF, /TNAME) EQ 'STRUCT' THEN BEGIN
			win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
			arr_Cp = ParametricCp_to_CurveCp(ST_AIF=*(*Pst).parametric_AIF, ARR_TIME=*(*Pst).results.arr_time, DOSE=(*Pst).par.dose, INJECTION_FRAME=(*Pst).par.nframe_injection)
			ok = Plot_profile_dynMRI(arr_Cp, ARR_TIME=*(*Pst).results.arr_time, WIN_NUMBER=win_number, METHOD=3)
		ENDIF
	END
	ELSE : BEGIN
		IF N_ELEMENTS(*(*Pst).results.arr_Signal_AIF) GT 1 THEN BEGIN
			arr_data = 	*(*Pst).results.arr_Signal_AIF
			arr_time =  *(*Pst).results.arr_time_AIF
			win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
			ok = Plot_profile_dynMRI(arr_data, ARR_TIME=arr_time, WIN_NUMBER=win_number, METHOD=2)
		ENDIF
	END
	ENDCASE

END
;-------------------------------------------------------------------------------------------------------
'PLOT_RR' : BEGIN
	CASE (*Pst).opt_RRcurve OF
	0 :
	ELSE : BEGIN
		IF N_ELEMENTS(*(*Pst).results.arr_Signal_RR) GT 1 THEN BEGIN
			arr_data = 	*(*Pst).results.arr_Signal_RR
			arr_time =  *(*Pst).results.arr_time_RR
			win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
			ok = Plot_profile_dynMRI(arr_data, ARR_TIME=arr_time, WIN_NUMBER=win_number, METHOD=2)
		ENDIF
	END
	ENDCASE

END
;-------------------------------------------------------------------------------------------------------
'BVIEW_RCE' : BEGIN

	window_mapped = WIDGET_INFO((*Pst).wd.wbases_images[1], /MAP)
	window_mapped_now = window_mapped EQ 0

	IF window_mapped_now EQ 1 THEN BEGIN
		ok = fwdcemri_Draw_RCEImage(Pst)
		sizeb1 = WIDGET_INFO((*Pst).wd.wbases_images[0], /GEOMETRY)
		WIDGET_CONTROL, (*Pst).wd.wviewbar_rce, SET_BUTTON=1
		WIDGET_CONTROL, (*Pst).wd.wbases_images[1], MAP=1, XSIZE=sizeb1.xsize, YSIZE=sizeb1.ysize
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*Pst).wd.wbases_images[1], MAP=0, XSIZE=1, YSIZE=1
		WIDGET_CONTROL, (*Pst).wd.wviewbar_rce, SET_BUTTON=0
	ENDELSE

END
;-------------------------------------------------------------------------------------------------------
'RCE_CALCULATION_TYPE' : BEGIN

IF ((*Pst).opt_typeRCE EQ 1) AND (ev.id EQ (*Pst).wd.wbttn_RCE_type1) THEN RETURN ; nothing.. (press already button set)
IF ((*Pst).opt_typeRCE EQ 2) AND (ev.id EQ (*Pst).wd.wbttn_RCE_type2) THEN RETURN ; nothing.. (press already button set)

IF ((*Pst).opt_typeRCE NE 1) AND (ev.id EQ (*Pst).wd.wbttn_RCE_type1) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbttn_RCE_type1, SET_BUTTON=1
	WIDGET_CONTROL, (*Pst).wd.wbttn_RCE_type2, SET_BUTTON=0
	(*Pst).opt_typeRCE = 1
ENDIF
IF ((*Pst).opt_typeRCE NE 2) AND (ev.id EQ (*Pst).wd.wbttn_RCE_type2) THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbttn_RCE_type1, SET_BUTTON=0
	WIDGET_CONTROL, (*Pst).wd.wbttn_RCE_type2, SET_BUTTON=1
	(*Pst).opt_typeRCE = 2
ENDIF

ok = fwdcemri_Draw_RCEImage(Pst, /CALCULATE)

END
;-------------------------------------------------------------------------------------------------------
'PLOTSOPENED_TYPE' : BEGIN

	(*Pst).opt_plotsOpened = fwdcemri_PlotsOpenedType(ev.id, (*Pst).wd.wbttn_plotsOpened, (*Pst).opt_plotsOpened)

END
;-------------------------------------------------------------------------------------------------------
'CHANGE_PALETTE' : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
		'PALETTE1' : BEGIN
			(*Pst).changing_palette = 1
			 (*Pst).vi1.opalette   -> Getproperty, RED=r,GREEN=G, BLUE=b
			 TVLCT, r,g,b
		END
		'PALETTE2' : BEGIN
			(*Pst).changing_palette = 2
			 (*Pst).vi2.opalette   -> Getproperty, RED=r,GREEN=G, BLUE=b
			 TVLCT, r,g,b
		END
	END
    XLOADCT, GROUP=(*Pst).wd.wbase_main, /MODAL, UPDATECALLBACK='fwdcemri_Change_palette', UPDATECBDATA=*Pst, /USE_CURRENT
END
;-------------------------------------------------------------------------------------------------------
'BACKGROUND_IMAGE_SIGNAL' : BEGIN

	WIDGET_CONTROL, (*Pst).wd.wDraw1, GET_VALUE = window_draw
	CASE (*Pst).opt.background_MRI OF
	1 : BEGIN ; first color of palette
		(*Pst).vi1.opalette->getproperty, RED_VALUES=r, GREEN_VALUES=g, BLUE_VALUES=b
		(*Pst).vi1.oview->setproperty, COLOR=[r[0],g[0],b[0]]
		(*Pst).vi1.ohcolorbar->setproperty, COLOR= (*Pst).opt.colorbar_labelcolor
		(*Pst).opt.background_MRI = 2
	END
	2 : BEGIN ; white color (and numbers of palette in black)
		(*Pst).vi1.oview->setproperty, COLOR=[255,255,255L]
		(*Pst).vi1.ohcolorbar->setproperty, COLOR= [0,0,0d]
		(*Pst).opt.background_MRI = 3
	END
	3 : BEGIN ; default color
		(*Pst).vi1.oview->setproperty, COLOR=(*Pst).opt.background_color_mri
		(*Pst).vi1.ohcolorbar->setproperty, COLOR= (*Pst).opt.colorbar_labelcolor
		(*Pst).opt.background_MRI = 1
	END
	ENDCASE

    window_draw -> Draw, (*Pst).vi1.oView
END
;-------------------------------------------------------------------------------------------------------
'BACKGROUND_IMAGE_RCE' : BEGIN

	WIDGET_CONTROL, (*Pst).wd.wDraw2, GET_VALUE = window_draw
	CASE (*Pst).opt.background_RCE OF
	1 : BEGIN ; first color of palette
		(*Pst).vi2.opalette->getproperty, RED_VALUES=r, GREEN_VALUES=g, BLUE_VALUES=b
		(*Pst).vi2.oview->setproperty, COLOR=[r[0],g[0],b[0]]
		(*Pst).vi2.ohcolorbar->setproperty, COLOR= (*Pst).opt.colorbar_labelcolor
		(*Pst).opt.background_RCE = 2
	END
	2 : BEGIN ; white color (and numbers of palette in black)
		(*Pst).vi2.oview->setproperty, COLOR=[255,255,255L]
		(*Pst).vi2.ohcolorbar->setproperty, COLOR= [0,0,0d]
		(*Pst).opt.background_RCE = 3
	END
	3 : BEGIN ; default color
		(*Pst).vi2.oview->setproperty, COLOR=(*Pst).opt.background_color_RCE
		(*Pst).vi2.ohcolorbar->setproperty, COLOR= (*Pst).opt.colorbar_labelcolor
		(*Pst).opt.background_RCE = 1
	END
	ENDCASE
	window_draw -> Draw, (*Pst).vi2.oView
END
;-------------------------------------------------------------------------------------------------------
'COLORBAR_FONT' : BEGIN

	uname= WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'INCREASE': val =  1
	'DECREASE': val = -1
	ENDCASE

	(*Pst).vi1.ohColorbar->Getproperty, FONTSIZE= fontsize
	(*Pst).vi1.ohColorbar->Setproperty, FONTSIZE= 1 > (fontsize+val)
	(*Pst).vi2.ohColorbar->Getproperty, FONTSIZE= fontsize
	(*Pst).vi2.ohColorbar->Setproperty, FONTSIZE= 1 > (fontsize+val)

	WIDGET_CONTROL, (*Pst).wd.wDraw1, GET_VALUE = window_draw1
	WIDGET_CONTROL, (*Pst).wd.wDraw2, GET_VALUE = window_draw2
	window_draw1 -> Draw, (*Pst).vi1.oView
	window_draw2 -> Draw, (*Pst).vi2.oView
END
;-------------------------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------------------------
'ZOOM_CHANGE' : BEGIN     ; Cambia el zoom

    uname= WIDGET_INFO(ev.id, /UNAME)
    WIDGET_CONTROL, (*Pst).wd.wbase_main, UPDATE=0
    CASE uname OF
    'ZOOM_MENOS' : BEGIN
    	ok = fwdcemri_Draw_Image(Pst, CONTROL_ZOOM=0.8d, NODRAW=1)
    	ok = fwdcemri_Draw_RCEImage(Pst, CONTROL_ZOOM=0.8d, NODRAW=1)
    END
    'ZOOM_MAS'   : BEGIN
    	ok = fwdcemri_Draw_Image(Pst, CONTROL_ZOOM=1.25d, NODRAW=1)
    	ok = fwdcemri_Draw_RCEImage(Pst, CONTROL_ZOOM=1.25d, NODRAW=1)
    END
    ENDCASE
    ok = fwdcemri_DrawInSameSize(Pst)
    ok = fwdcemri_REDraw(Pst)
    WIDGET_CONTROL, (*Pst).wd.wbase_main, UPDATE=1 ;

END
;-------------------------------------------------------------------------------------------------------
'INTERPOLATION': BEGIN
     (*Pst).flag_interp = WIDGET_INFO((*Pst).wd.wbttn_interp, /BUTTON_SET)
     ok=fwdcemri_Draw_Image(Pst, NODRAW=0)
     ok=fwdcemri_Draw_RCEImage(Pst, NODRAW=0)
     ;ok = fwdcemri_REDraw(Pst)
END
;-------------------------------------------------------------------------------------------------------
'ROI_SELECTION': BEGIN
	(*Pst).flag_drawROI = WIDGET_INFO((*Pst).wd.wbttn_drawROI, /BUTTON_SET)
	CASE (*Pst).flag_drawROI OF
	0 : ok = fwdcemri_activate_bases(Pst,  /NO_DRAW_ROI)
    1 : ok = fwdcemri_activate_bases(Pst, /DRAW_ROI)
     ENDCASE

    ok = fwdcemri_Redraw(Pst)
	(*Pst).flag_noresetroi = 0

END
;-------------------------------------------------------------------------------------------------------
'NEW_ROI': BEGIN

	ok = fwdcemri_Activate_Bases(Pst, /NEW_ROI)

	(*Pst).flag_analysis      = 0
	(*Pst).flag_analysis_done = 0

	PRINT, (*Pst).flag_noresetroi
	IF (*Pst).flag_noresetroi EQ 0 THEN BEGIN
		ok = fwdcemri_resetROI(Pst)
	ENDIF
	(*Pst).flag_noresetroi = 0

	(*Pst).vi1.ohcolorbar-> setproperty, PALETTE=(*Pst).vi1.opalette, RANGE=[0,255]
	;(*Pst).vi1.oimageROI -> setProperty, HIDE=1

	(*Pst).vi1.Opolyline_ROI ->SetProperty, HIDE=1
	(*Pst).vi2.Opolyline_ROI ->SetProperty, HIDE=1
	(*Pst).vi1.Opolyline_AUX ->SetProperty, HIDE=1
	(*Pst).vi2.Opolyline_AUX ->SetProperty, HIDE=1

	ok = fwdcemri_Redraw(Pst)

END
;-------------------------------------------------------------------------------------------------------
'ROI_TYPE' : BEGIN
	(*Pst).roi.opt_typeROI = WIDGET_INFO((*Pst).wd.wdropl_typeROI, /DROPLIST_SELECT)
	(*Pst).flag_drawROI    = WIDGET_INFO((*Pst).wd.wbttn_drawROI, /BUTTON_SET)

	ok = fwdcemri_resetROI(Pst)
	CASE (*Pst).roi.opt_typeROI OF
	0 : BEGIN ; box
		WIDGET_CONTROL, (*Pst).wd.wtext_ROI_size[0], SEND_EVENT={ID:0L, TOP:0L, HANDLER:0L}
		(*Pst).flag_analysis = 1 ; directly can be analized
	END
	1 : BEGIN ; whole image
		point_ini = fwdcemri_fullROI_set(Pst, DRAWTWO=1, ID_TEXT=[(*Pst).wd.wtext_info1,(*Pst).wd.wtext_info2])
		*(*Pst).roi.idx_points = IdxPointsFromBox(POINT_INI=[0,0l], SIZE_ROI=(*Pst).inf.size_im, DIMENSIONS=(*Pst).inf.size_im)
		(*Pst).flag_analysis = 1 ; directly can be analized
	END
	2 : BEGIN ; free
		(*Pst).flag_analysis = (*Pst).flag_freeROIclose
	END
	ENDCASE
	ok = fwdcemri_redraw(Pst)
	ok = fwdcemri_activate_bases(Pst, /CHANGE_ROI_TYPE)

END
;-------------------------------------------------------------------------------------------------------
'TEXT_ROI' : BEGIN

	uname= WIDGET_INFO(ev.id, /UNAME)
	WIDGET_CONTROL, (*Pst).wd.wtext_roi_px, GET_VALUE=str_point_x
	IF str_point_x EQ '' THEN RETURN
	WIDGET_CONTROL, (*Pst).wd.wtext_roi_py, GET_VALUE=str_point_y
	IF str_point_y EQ '' THEN RETURN

	ok = fwdcemri_ROI_sizeresolution(Pst) ; read size and resolution

	point_ini_x = 0 > FIX(str_point_x[0]) < ((*Pst).inf.size_im[0]-1)
	point_ini_y = 0 > FIX(str_point_y[0]) < ((*Pst).inf.size_im[1]-1)

	point_ini = get_point_YdownConvention([point_ini_x,point_ini_y],$
		DIMENSIONS=(*Pst).inf.size_im, SIZE_ROI=(*Pst).inf.size_roi)
	; ¡se guarda point_ini en IDL convention!

	point_ini[0] = 0 > point_ini[0] < ((*Pst).inf.size_im[0]-(*Pst).inf.size_roi[0])
	point_ini[1] = 0 > point_ini[1] < ((*Pst).inf.size_im[1]-(*Pst).inf.size_roi[1])

	(*Pst).inf.p_ini = point_ini

	; only with non free ROIs change idx_points
	CASE (*Pst).roi.opt_typeROI OF
	0 : BEGIN
		point_ini = fwdcemri_boxROI_set(Pst, POINT_INI=point_ini, DRAWTWO=1, $
			ID_TEXT=[(*Pst).wd.wtext_info1,(*Pst).wd.wtext_info2])

		IF N_ELEMENTS(point_ini) NE 2 THEN BEGIN ; error, hay que borrar el punto marcado
			CASE uname OF
				'POINT_ROI_X': WIDGET_CONTROL, (*Pst).wd.wtext_roi_px, SET_VALUE=''
				'POINT_ROI_Y': WIDGET_CONTROL, (*Pst).wd.wtext_roi_py, SET_VALUE=''
				ELSE:
			ENDCASE
		ENDIF
		*(*Pst).roi.idx_points = IdxPointsFromBox(POINT_INI=point_ini, SIZE_ROI=(*Pst).inf.size_roi, $
			DIMENSIONS=(*Pst).inf.size_im)
	END
	1 : BEGIN
		point_ini = fwdcemri_fullROI_set(Pst, DRAWTWO=1, ID_TEXT=[(*Pst).wd.wtext_info1,(*Pst).wd.wtext_info2])
		*(*Pst).roi.idx_points = IdxPointsFromBox(POINT_INI=[0,0l], SIZE_ROI=(*Pst).inf.size_im, $
			DIMENSIONS=(*Pst).inf.size_im)
	END
	2 :
	ENDCASE
	;PRINT, 'Resolucion:', (*Pst).roi.np_lr
	ok = fwdcemri_redraw(Pst)
	ok = fwdcemri_Activate_bases(Pst, /ALLOW_ANALYSIS)

END
;--------------------------------------------------------------------------------------------------
'OPEN_MASK_TXT' : BEGIN
	WIDGET_CONTROL, (*Pst).wd.wDraw1,  GET_VALUE = window_1
	window_1->getproperty, IMAGE_DATA=image_view

	IF SIZE(image_view, /N_DIMENSIONS) NE 3 THEN BEGIN
		ok = DIALOG_MESSAGE('Data must be loaded first', /INFO) &	RETURN
	ENDIF

	str_file = DIALOG_PICKFILE(PATH=(*Pst).path_result, TITLE='Open ROI txt file', DIALOG_PARENT=ev.top, $
		DEFAULT_EXTENSION='TXT')
	IF FILE_TEST(str_file) NE 1 THEN RETURN

	str_array = ReadFile_ASCII(str_file, REMOVE_BLANKS=1, ERROR_ID=error_id)
	IF error_id NE 0 THEN BEGIN
		ok = DIALOG_MESSAGE('Error loading txt file', /INFO) &	RETURN
	ENDIF
	im_mask   = Get_ROIfromTxTdata(str_array, DIMENSIONS=dimensions, RESOLUTION=resolution, ERROR_ID=error_id)
	IF error_id NE 0 THEN BEGIN
		ok = DIALOG_MESSAGE('Error in file format', /INFO) &	RETURN
	ENDIF
	IF (dimensions[0] NE (*Pst).inf.size_im[0]) OR (dimensions[1] NE (*Pst).inf.size_im[1]) THEN BEGIN
		ok = DIALOG_MESSAGE('ROI dimensions not equal to dimensions of current image', /INFO) &	RETURN
	ENDIF

	;window_1->show, 0
	(*Pst).vi1.oimage->getproperty, DIMENSIONS=dimensions2, SUB_RECT=sub_rect ; dimensions2 must be equal to dim

	img_color = BYTARR(dimensions[0],dimensions[1],3)
	mask_alphachannel = [[[img_color]],[[(im_mask EQ 0)*255B]]]
	(*Pst).vi1.oimageMASK -> setProperty, DATA = mask_alphachannel, DIMENSIONS=dimensions, GREYSCALE=0,$
    	SUB_RECT = sub_rect, HIDE=0

    window_1-> Draw, (*Pst).vi1.oView
    ok = fwdcemri_Activate_bases(Pst, /OPEN_MASK)

END
;--------------------------------------------------------------------------------------------------
'CLOSE_MASK_TXT' : BEGIN
	WIDGET_CONTROL, (*Pst).wd.wDraw1,  GET_VALUE = window_1
	(*Pst).vi1.oimageMASK -> setProperty, HIDE=1
	window_1-> Draw, (*Pst).vi1.oView
	ok = fwdcemri_Activate_bases(Pst, /FREE_MASK)
END
;--------------------------------------------------------------------------------------------------
'EXPORT_3DSTACK' : BEGIN ; provisional para grabar por separados los volumenes dinamicos

;data = (*(*Pst).data)

size_datam = [(*Pst).inf.size_im[0],(*Pst).inf.size_im[1],(*Pst).inf.nslices,(*Pst).inf.nframes]

str_fileseq = 'F:\Experiments_DCE-MRI\dMRI_milena110524\7\pdata\1\2dseq'

data = Read_Binary_BrukerData(str_fileseq, ROTATION=0,$
    	DIMENSIONS=size_data, TYPE_DATA=2) ; integer

length_numbers = STRLEN(STRTRIM((*Pst).inf.nframes,2))

;dim = SIZE(data_int, /DIMENSIONS)
;strinfo_size   = STRTRIM(dim[0],2) + 'x' +  STRTRIM(dim[1],2) +  'x' + STRTRIM(dim[2],2)
;str_name = 'Vol_float_' + strinfo_size

str_filebase = DIALOG_PICkFILE(PATH=(*Pst).path_result, TITLE='save 3D stack', DIALOG_PARENT=ev.top, $
		DEFAULT_EXTENSION='', /WRITE)

;'vol_128x128x6_int_'

IF str_filebase EQ '' THEN RETURN

FOR i=0l, (*Pst).inf.nframes-1 DO BEGIN
	str_nframe = Translate_number2string(NUMBER=i+1, LENGTH=length_numbers)
	str_file = str_filebase + '_' + str_nframe + '.bin'
	vol = REFORM(data[*,*,*,i])
	ok = write_raw(FILE=str_file, DATA=vol)
ENDFOR

END
;--------------------------------------------------------------------------------------------------
'EXPORT_WINDOW' : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'EXPORT_WIN1' : BEGIN
		WIDGET_CONTROL, (*Pst).wd.wDraw1,  GET_VALUE = window_1
		;(*Pst).vi1.owindow->getproperty, IMAGE_DATA=image_view
	END
	'EXPORT_WIN2' : BEGIN
		WIDGET_CONTROL, (*Pst).wd.wDraw2,  GET_VALUE = window_1
		;(*Pst).vi2.owindow->getproperty, IMAGE_DATA=image_view
	END
	ENDCASE
	window_1->getproperty, IMAGE_DATA=image_view
	path_result = (*Pst).path_result
	im = fwdcemri_exportImage(image_view, OPTION_NAME=0, PATH_RESULT=path_result, INFO= *(*Pst).inf.info)
	(*Pst).path_result = path_result

END
;-------------------------------------------------------------------------------------------------------
'ANALYZE_ROI' : BEGIN

	modelnumber = WIDGET_INFO((*Pst).wd.wdropls_model, /DROPLIST_SELECT)

	st_par = (*Pst).par
	ok= fwdcemri_Parameters_read (ST_WIDGETS=(*Pst).wdp, ST_PARAM=st_par, ERROR_ID=error_id, ERROR_STR=error_str)
	ok= fwdcemri_Parameters_write(ST_WIDGETS=(*Pst).wdp, ST_PARAM=st_par)
	IF error_id NE 0 THEN BEGIN
		ok = DIALOG_MESSAGE(error_str, /INFO) &	RETURN ; no es una estructura
	ENDIF
	(*Pst).par = st_par

	xsize     = (*Pst).inf.size_im[0]
	ysize     = (*Pst).inf.size_im[1]
	nframes   = (*Pst).inf.nframes
	n_z       = (*Pst).inf.current_slice

	roi_st = UpdateROIdata(IDX_POINTS=*(*Pst).roi.idx_points, ST_ROI=(*Pst).roi, IMAGESIZE=(*Pst).inf.size_im, SLICEZ=n_z, ERROR_ID=error_id, ERROR_STR=error_str)
	IF error_id NE 0 THEN BEGIN
		ok = DIALOG_MESSAGE(error_str, /INFO) & RETURN & ENDIF
	(*Pst).roi = roi_st
	;-------------------------------------------------------------------------------------------------------
	pos_widget_child = locate_FirstNoValid_Widget((*Pst).wd.wAuxMainBases)
	IF pos_widget_child EQ -1 THEN BEGIN
		ok = DIALOG_MESSAGE('Maximun 10 opened simultaneous results', /INFO)
		RETURN
	ENDIF
	;-------------------------------------------------------------------------------------------------------

    *(*Pst).results.arr_time = INDGEN(st_par.nframes)*st_par.frame_period
	*(*Pst).modelinfo = DCEMRI_ConfigModelParameters(modelnumber, ERROR_ID=error_id)
	;-------------------------------------------------------------------------------------------------------
	IF error_id NE 0 THEN BEGIN
		ok = DIALOG_MESSAGE('Internal error #001, please contact software provider', /INFO) & RETURN & ENDIF
	;-------------------------------------------------------------------------------------------------------

	IF (*Pst).opt_AIFcurve EQ 1 THEN BEGIN
		IF SIZE(*(*Pst).parametric_aif, /TNAME) EQ 'STRUCT' THEN BEGIN ; it's biexponential
			*(*Pst).parametric_aif= fwdcemri_ParametricAIF_read( ST_WIDGETS=(*Pst).wdp, ST_AIF=*(*Pst).parametric_aif, ERROR_ID=error_id)
		ENDIF
	ENDIF

	out_results = fwdcemri_Perform_ModelAnalysis(ST_RESULTS=(*Pst).results, DATA=*(*Pst).data, ST_AIF=*(*Pst).parametric_AIF, ST_MODELINFO=*(*Pst).modelinfo,$
		ST_ROI=(*Pst).roi, ST_PAR=(*Pst).par, TYPE_AIF=(*Pst).opt_AIFcurve, ERROR_ID=error_id, T1MAP=*(*Pst).t1map)

	IF error_id EQ 0 THEN BEGIN
		(*Pst).results = out_results
		Interface_DCEMRIresult, pst, TITLE=(*(*Pst).modelinfo).str_model, WIDGET_ID=widget_id_mainbase
		(*Pst).wd.wAuxMainBases[pos_widget_child] =widget_id_mainbase
	ENDIF


END
;-------------------------------------------------------------------------------------------------------
'DRAW' : BEGIN     ; Drawing events
	ok = fwdcemri_DrawingEvents(ev, pst)
END
;-------------------------------------------------------------------------------------------------------
'BABOUT' : BEGIN
	ok = Interface_DCEMRI_About(GROUP_LEADER=(*Pst).wd.wbase_main)
	RETURN
END
;-------------------------------------------------------------------------------------------------------
'BQUIT' : BEGIN       ; Finaliza el programa
	ok = DIALOG_MESSAGE('Quit?', /QUESTION)
	IF STRUPCASE(ok) EQ 'NO' THEN RETURN
    Interface_DCEMRI_destroy, ev.top
    HEAP_GC &  RETURN
END ;BQUIT
;-------------------------------------------------------------------------------------------------------
ELSE: PRINT, 'Identificador inválido'
ENDCASE


END ;Interface_DCEMRI_event

;**************************************************************************************************
;**************************************************************************************************

PRO  Interface_DCEMRI, data, TITLE=title, FILENAME=filename, FOV=fov, IMAGE_BITMAP=image_bitmap, $
	WIDGET_ID = widget_id

COMPILE_OPT DEFINT32, STRICTARR ; , HIDDEN

opt_catcherror = 0
;*****************************************
IF LMGR(/RUNTIME) OR LMGR(/VM) OR (opt_catcherror EQ 1) THEN BEGIN
	Catch, theError
	IF theError NE 0 THEN BEGIN
		Catch, /Cancel
		void = DIALOG_MESSAGE(!ERROR_STATE.MSG, /ERROR)
		;void = ERROR_MESSAGE()
		RETURN
	ENDIF
ENDIF
;*****************************************

title_base_main = 'DCE@urLAB'
str_file_config = 'DCEurLAB_config.txt'

n_instance =  XRegistered('Interface_DCEMRI')
IF n_instance GT 0 THEN title_base_main+= + ' (' + STRTRIM(FIX(n_instance[0]+1),2) + ')'

IF FILE_TEST(str_file_config,/REGULAR) NE 1 THEN $
	str_file_config = DIALOG_PICKFILE(TITLE='Choose Config file (DCEurLAB_config.txt) file is missing in main directory')

err_read = 0l
;----------------------------------------------------------
IF FILE_TEST(str_file_config ) EQ 1 THEN BEGIN
	array_in  = ReadFile_ASCII(str_file_config, REMOVE_BLANKS=1, ERROR_ID=error_id, MAX_LINES=1e5, MAX_BYTES=1e6)
	IF error_id NE  0 THEN BEGIN
		ok = DIALOG_MESSAGE('Error reading config file (1), load default values and continue?', /QUESTION)
		IF STRUPCASE(ok) EQ 'NO' THEN RETURN
		err_read = 1l
	ENDIF ELSE BEGIN
		array_in = Remove_Comments(array_in, STR_COMMENTS='#')
		st_param = DCEMRI_ReadConfig(ARRAY_IN=array_in, ARRAY_OUT=array_out, $
			ERROR_ID=error_id, ERROR_STR=error_str, ST_AUX=stp, ST_AIF=st_aif)
		IF error_id NE 0 THEN BEGIN
			ok = DIALOG_MESSAGE([['Error reading config file (2) '],[error_str],[''],['load default values and continue?']], /QUESTION)
			IF STRUPCASE(ok) EQ 'NO' THEN RETURN
			err_read = 2
		ENDIF
	ENDELSE
ENDIF ELSE BEGIN
	ok = DIALOG_MESSAGE([['Missing configuration file'],['load default values and continue?']], /QUESTION)
	IF STRUPCASE(ok) EQ 'NO' THEN RETURN
	err_read = 3
ENDELSE
;----------------------------------------------------------
IF err_read NE 0 THEN BEGIN

	ok = DIALOG_MESSAGE('Default set of values will be used', /INFO)
	; default parameters
	st_param = DCEMRI_InitialParameters(ST_AIF=st_aif)
	stp = {sizeview:360l,$
		flag_interpolation : 0,$
		nframes  :st_param.nframes,$
		nslices  :4l,$
		size_roi : [40,40l],$
		np_lr    : 1l, $
		type_roi : 2l,$
		size_im  : [256,256l],$
		path_data:'',path_t1map:'', path_result:'', path_aif:''}

		CD, CURRENT=path_unit
		stp.path_Data    =path_unit
		stp.path_t1map  = path_unit
		stp.path_result = path_unit
		stp.path_aif    = path_unit
ENDIF

;----------------------------------------------------------------------
st_inf  = { $
    size_im :  stp.size_im,  $
	nslices :  stp.nslices,  $
	nframes :  stp.nframes,  $ ;
	current_slice         : 0l,$
	current_frame         : 0l,$
	slope   :  PTR_NEW(/ALLOCATE_HEAP), $ ; In Bruker format, value to multiply data read in disk
	info    :  PTR_NEW(/ALLOCATE_HEAP), $ ; Additional information about the image (name, directory, etc...)
	file    : '',    $
	name    : '',	 $
	path    : '',	 $
	size_t2 : [200l,200l], $
	nslices_T2 :   3L,  $
	slope_T2   : 1.0,   $
	file_T2    : '',    $
	p_ini      :  [-1,-1], $    ;   initial point, IDL convention (left-down corner)
	size_roi   :  stp.size_roi $;  In case of free ROI, size of the rectangle covering the ROI

}

*st_inf.slope = 1.0 ; al principio el "slope" es 1

DEVICE, Get_Screen_Size = screenSize
view_size_limits_X = [200, screenSize[0]/3]   ; mayor y menor en X
view_size_limits_y = [200, screenSize[1]-300] ; mayor y menor en Y
view_scroll_limit = 2500l ; absolute limit for image even with scrolling
;---------------------------------------------------------------------
st_draw_options = { $
	size_view_ini    : stp.sizeview, $
	draw_size        :INTARR(2),$
	draw_scroll_size :INTARR(2),$
	background_RCE : 1, $
	background_MRI : 1, $
	background_color_RCE  : [100,100,100], $
	background_color_mri  : [100,100,100], $
	; values to configure the view window.
	view_size_limits_x      : view_size_limits_x,$
	view_size_limits_y      : view_size_limits_y,$
	view_scroll_limit       : view_scroll_limit, $
	val_frame1            : 0.125/4.0, $ ; Frame around the image
	scale_relation_xy     : 0d,$    ; Very important. Its an additional compensation for the colorbar
	position_colorbar_ini : FLTARR(4),$
	colorbar_relativesize : 0.30,$  ; Additional reserved space for colorbar, relative to image size
	colorbar_labelcolor   : [255,255,255],$
	colorbar_fontsize     : 14, $
	info : 0l  $
}
st_draw_options.scale_relation_xy = (1/(st_draw_options.colorbar_relativesize+1))
st_draw_options.draw_size = [st_draw_options.size_view_ini, st_draw_options.size_view_ini*(1+st_draw_options.colorbar_relativesize[0])] ; take into account the colorbar (in scale relation_xy)
st_draw_options.position_colorbar_ini = [0.1, 0.08, 0.9, 0.98]
st_draw_options.position_colorbar_ini[1]+=st_draw_options.scale_relation_xy
st_draw_options.draw_scroll_size = st_draw_options.draw_size
;---------------------------------------------------------------------

IF N_ELEMENTS(title) NE 0 THEN title_base_main+=title[0]
IF N_ELEMENTS(fov) NE 3 THEN size_fov = DBLARR(3) ELSE size_fov = fov[0:2]*1d

idx_paletteRCE_ini = 33 ; bluered colormap
idx_paletteROI_ini = 33 ; bluered colormap

size_slider_1 = stp.sizeview*0.5
size_slider_2 = stp.sizeview*0.25
size_slider_3 = 80

initial_zoom = 1

base_main  =WIDGET_BASE(/COLUMN, TITLE=title_base_main, MBAR=barBase, TLB_FRAME_ATTR=1, MAP=0,$
	BITMAP=image_bitmap, KILL_NOTIFY='Interface_DCEMRI_killnotify')

widget_id = base_main

base_smain =WIDGET_BASE(base_main, /COLUMN, FRAME=0)
base_a     =WIDGET_BASE(base_smain, /ROW,  FRAME=0,   /ALIGN_RIGHT, /BASE_ALIGN_RIGHT)
	base_a0  =WIDGET_BASE(base_a, /ROW,  FRAME=0,   /ALIGN_RIGHT, /BASE_ALIGN_RIGHT)
base_b     =WIDGET_BASE(base_smain, /ROW, FRAME=0)
	base_b1    =WIDGET_BASE(base_b, /COLUMN,  FRAME=0,  /ALIGN_LEFT, /BASE_ALIGN_CENTER, YPAD=0, XPAD=0)
	base_b1a=WIDGET_BASE(base_b1, /COLUMN)
		base_rell = WIDGET_BASE(base_b1, YSIZE=5)
	base_b1b=WIDGET_BASE(base_b1, /COLUMN, /ALIGN_CENTER, /BASE_ALIGN_CENTER, /FRAME, YPAD=5)
		base_b1b1= WIDGET_BASE(base_b1b, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER, FRAME=0)
			base_b1b1a = WIDGET_BASE(base_b1b1, /ROW,  /ALIGN_CENTER, /BASE_ALIGN_CENTER, /NONEXCLUSIVE)
			base_b1b1b = WIDGET_BASE(base_b1b1, /ROW,  /ALIGN_CENTER, /BASE_ALIGN_CENTER)
			base_b1b1c = WIDGET_BASE(base_b1b1, /ROW,  /ALIGN_CENTER, /BASE_ALIGN_CENTER)
		base_b1b2= WIDGET_BASE(base_b1b, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER, FRAME=0)
			base_b1b2a = WIDGET_BASE(base_b1b2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
			base_b1b2b = WIDGET_BASE(base_b1b2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_b1b3= WIDGET_BASE(base_b1b, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
			base_b1b3a = WIDGET_BASE(base_b1b3, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
			base_b1b3b = WIDGET_BASE(base_b1b3, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
	base_rell = WIDGET_BASE(base_b1, YSIZE=5)
	base_b1c=WIDGET_BASE(base_b1, /COLUMN, /ALIGN_CENTER, /BASE_ALIGN_CENTER, FRAME=1)
		base_b1c1 = WIDGET_BASE(base_b1c, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
		base_b1c2 = WIDGET_BASE(base_b1c, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_b1c3 = WIDGET_BASE(base_b1c, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_b1c4 = WIDGET_BASE(base_b1c, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
	base_rell = WIDGET_BASE(base_b1, YSIZE=5)
	base_b1d=WIDGET_BASE(base_b1, /COLUMN, /ALIGN_CENTER, /BASE_ALIGN_CENTER, FRAME=1)
		base_b1d1 = WIDGET_BASE(base_b1d, /COLUMN, /ALIGN_CENTER, /BASE_ALIGN_CENTER)

base_rell = WIDGET_BASE(base_b, /COLUMN, FRAME=0, XSIZE=5)
base_b2 =WIDGET_BASE(base_b, /COLUMN, FRAME=1, /ALIGN_CENTER, YPAD=5, XPAD=5)
	base_b2a =WIDGET_BASE(base_b2, /COLUMN, FRAME=0, /ALIGN_CENTER, YPAD=0, XPAD=0)
		base_rell=WIDGET_BASE(base_b2a, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER, YSIZE=5)
 	base_b2b =WIDGET_BASE(base_b2, /COLUMN, FRAME=0, /ALIGN_CENTER)
	base_b2c =WIDGET_BASE(base_b2, /COLUMN, FRAME=0, /ALIGN_CENTER)
base_rell = WIDGET_BASE(base_b, /COLUMN, FRAME=0, XSIZE=5)
base_b3 =WIDGET_BASE(base_b, /COLUMN, FRAME=1, /ALIGN_CENTER, YPAD=5, XPAD=5)
	base_b3a =WIDGET_BASE(base_b3, /COLUMN, FRAME=0, /ALIGN_CENTER, YPAD=0, XPAD=0)
		base_rell=WIDGET_BASE(base_b3a, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER, YSIZE=5)
	base_b3b =WIDGET_BASE(base_b3, /COLUMN, FRAME=0, /ALIGN_CENTER)
	base_b3c =WIDGET_BASE(base_b3, /COLUMN, FRAME=0, /ALIGN_CENTER)

Filebar = WIDGET_BUTTON(barBase, VALUE='File', UVALUE='FILE', /MENU)
    Openbar_t1  = WIDGET_BUTTON(Filebar, VALUE='Open DCE-MRI', UVALUE='BOPEN_MRI', UNAME='OPEN_MRI', SENSITIVE=1, /MENU)
    Openbar_t1a = WIDGET_BUTTON(Openbar_t1, VALUE='Raw data', UVALUE='BOPEN_MRI_DT1', UNAME='OPEN_MRI_RAW', SENSITIVE=1)
    Openbar_t1b = WIDGET_BUTTON(Openbar_t1, VALUE='Bruker format', UVALUE='BOPEN_MRI_DT1', UNAME='OPEN_MRI_BRUKER', SENSITIVE=1)
    Openbar_t1c = WIDGET_BUTTON(Openbar_t1, VALUE='DICOM format', UVALUE='BOPEN_MRI_DT1',  UNAME='OPEN_MRI_DICOM', SENSITIVE=1)
    Openbar_t1free = WIDGET_BUTTON(Filebar, VALUE='Close DCE-MRI', UVALUE='BFREE_DATA', $
       UNAME='FREE_MRI_T1', SENSITIVE=0)
    Openbar_t1map      = WIDGET_BUTTON(Filebar, VALUE='Open T10 map', UVALUE='BOPEN_T1MAP',   UNAME='OPEN_T1MAP', SENSITIVE=1, /SEPARATOR)
    Openbar_t1mapfree  = WIDGET_BUTTON(Filebar, VALUE='Close T10 map', UVALUE='BFREE_DATA', UNAME='FREE_T1MAP', SENSITIVE=0)
    ;Openbar_mask      = WIDGET_BUTTON(Filebar, VALUE='Load mask', UVALUE='OPEN_MASK_TXT', SENSITIVE=0, /SEPARATOR)
    ;Openbar_maskfree  = WIDGET_BUTTON(Filebar, VALUE='Close mask', UVALUE='CLOSE_MASK_TXT', SENSITIVE=0)
	Openbar_AIF        = WIDGET_BUTTON(Filebar, VALUE='Open AIF', /MENU, /SEPARATOR, SENSITIVE=1)
        Openbar_AIF_curve = WIDGET_BUTTON(Openbar_AIF, VALUE='Curve file', UVALUE='OPEN_CURVE_TXT', UNAME='AIF_CURVE')
        Openbar_AIF_param = WIDGET_BUTTON(Openbar_AIF, VALUE='Parametric function file', UVALUE='OPEN_CURVE_PARAMETRIC', UNAME='AIF_PARAMETRIC')
	;Openbar_AIFfree   = WIDGET_BUTTON(Filebar, VALUE='Close AIF', UVALUE='CLOSE_AIF_TXT', SENSITIVE=0)
	Openbar_RR        = WIDGET_BUTTON(Filebar, VALUE='Open RR curve file', UVALUE='OPEN_CURVE_TXT', UNAME='RR_CURVE', SENSITIVE=1)
    ;Openbar_RRfree    = WIDGET_BUTTON(Filebar, VALUE='Close RR', UVALUE='CLOSE_RR_TXT', SENSITIVE=0)

    Quitbar = WIDGET_BUTTON(Filebar, VALUE='Quit', UVALUE='BQUIT', /SEPARATOR)

ExportBar = WIDGET_BUTTON(barBase, VALUE='Export\Import', /MENU, SENSITIVE=0)
	ExportWindow  = WIDGET_BUTTON(exportbar,     VALUE='Export images',  /MENU)
	exportbar_w1  = WIDGET_BUTTON(ExportWindow, VALUE='left Window',  UVALUE='EXPORT_WINDOW', UNAME='EXPORT_WIN1')
	exportbar_w2  = WIDGET_BUTTON(ExportWindow, VALUE='right Window', UVALUE='EXPORT_WINDOW', UNAME='EXPORT_WIN2')
	;Export3dstack  = WIDGET_BUTTON(exportbar,     VALUE='Export 3D image stack',  /SEPARATOR, UVALUE='EXPORT_3DSTACK')
	bttn_ExportROIkinetics = WIDGET_BUTTON(exportbar, VALUE='Export ROI kinetic curve', /SEPARATOR, /MENU)
		bttn_ExportROIkinetics_RCE = WIDGET_BUTTON(bttn_ExportROIkinetics, VALUE='RCE values (%)', UVALUE=['EXPORT_ASCII','COLUMN_FORMAT'], UNAME='RCE_VALUES')
		bttn_ExportROIkinetics_ABS = WIDGET_BUTTON(bttn_ExportROIkinetics, VALUE='MRI values (signal)', UVALUE=['EXPORT_ASCII','COLUMN_FORMAT'], UNAME='ABS_VALUES')
	bttn_ExportRoi = WIDGET_BUTTON(exportbar, VALUE='Export ROI', /SEPARATOR,	UVALUE='EXPORT_ROI', UNAME='EXPORT_ROI_SAV')
	bttn_importRoi = WIDGET_BUTTON(exportbar, VALUE='Import ROI', UVALUE='IMPORT_ROI', UNAME='IMPORT_ROI_SAV')

Viewbar = WIDGET_BUTTON(barBase, VALUE='View', UVALUE='VIEW', /MENU, SENSITIVE=0)
    viewbar_rce    = WIDGET_BUTTON(viewbar, VALUE='View RCE image', UVALUE='BVIEW_RCE', UNAME='VIEW_RCE',SENSITIVE=1, /CHECKED_MENU)
    viewbar_t1map  = WIDGET_BUTTON(viewbar, VALUE='View T10 map', UVALUE='BVIEW_T1MAP',  UNAME='VIEW_T1MAP',  SENSITIVE=0)
    viewbar_data1  = WIDGET_BUTTON(viewbar, VALUE='View dynamic frames', UVALUE='BVIEW_DYN_MRI', UNAME='VIEW_DYN_MRI',SENSITIVE=0, /SEPARATOR)
    viewbar_data2  = WIDGET_BUTTON(viewbar, VALUE='View cine', UVALUE='BVIEW_DYN_MRI', UNAME='VIEW_CINE',SENSITIVE=0)
    viewbar_ROIAIF = WIDGET_BUTTON(viewbar, VALUE='Plot ROI kinetic curve', UVALUE='ROI_KINETICS', UNAME='TYPE_PLOT_1',  SENSITIVE=0, /SEPARATOR)
    viewbar_ROIkinetics = WIDGET_BUTTON(viewbar, VALUE='Plot ROI kinetic curve (inc. SD)', UVALUE='ROI_KINETICS', UNAME='TYPE_PLOT_3',  SENSITIVE=0)
    viewbar_deleteWindows = WIDGET_BUTTON(viewbar, VALUE='Close secondary windows', UVALUE='CLOSE_WINDOWS', /SEPARATOR)

bttns_viewbarROIs = [viewbar_ROIkinetics,viewbar_ROIAIF]

Optionsbar  = WIDGET_BUTTON(barBase, VALUE='Options', UVALUE='OPTIONS', /MENU, SENSITIVE=0)

	barmenu_datainfo    = WIDGET_BUTTON(Optionsbar, VALUE='Image info', UVALUE='IMAGE_INFO')
	barmenu_limitframes = WIDGET_BUTTON(Optionsbar, VALUE='Select dynamic range', UVALUE='DYNAMIC_RANGE', UNAME='DYNAMIC_RANGE')
	barmenu_palettes = WIDGET_BUTTON(Optionsbar, VALUE='Color palette', /MENU, /SEPARATOR)
	    bar_palette1 = WIDGET_BUTTON(barmenu_palettes, VALUE='MRI', UVALUE='CHANGE_PALETTE', UNAME='PALETTE1')
	    bar_palette2 = WIDGET_BUTTON(barmenu_palettes, VALUE='RCE', UVALUE='CHANGE_PALETTE', UNAME='PALETTE2')
	bar_background_signal = WIDGET_BUTTON(Optionsbar, VALUE='Change background (MRI)', 	UVALUE='BACKGROUND_IMAGE_SIGNAL', /SEPARATOR)
	bar_background_rce    = WIDGET_BUTTON(Optionsbar, VALUE='Change background (RCE)', 	UVALUE='BACKGROUND_IMAGE_RCE')
	bar_colorbar_font     = WIDGET_BUTTON(Optionsbar, VALUE='Change colorbar font', 	UVALUE='COLORBAR_FONT', /MENU)
		bar_colorbar_mas  = WIDGET_BUTTON(bar_colorbar_font, VALUE='Increase size', 	UVALUE='COLORBAR_FONT', UNAME='INCREASE')
		bar_colorbar_menos= WIDGET_BUTTON(bar_colorbar_font, VALUE='Decrease size', 	UVALUE='COLORBAR_FONT', UNAME='DECREASE')
	bar_initialparam  = WIDGET_BUTTON(Optionsbar, VALUE='Initial parameters',   UVALUE='INITIAL_PARAMETERS', /SEPARATOR)

	bar_RCE_type = WIDGET_BUTTON(Optionsbar, VALUE='RCE calculation', /MENU, /SEPARATOR)

	bttn_RCE_type1 = WIDGET_BUTTON(bar_RCE_type, VALUE='Relative to values before CA injection', /CHECKED_MENU, UVALUE='RCE_CALCULATION_TYPE', UNAME='RCE_TYPE_1')
	bttn_RCE_type2 = WIDGET_BUTTON(bar_RCE_type, VALUE='Absolute (max/min)',  /CHECKED_MENU,     UVALUE='RCE_CALCULATION_TYPE', UNAME='RCE_TYPE_2')
	bar_PlotsOpened_type = WIDGET_BUTTON(Optionsbar, VALUE='View plots', /MENU, /SEPARATOR)
		bttn_PlotsOpened_type1 =  WIDGET_BUTTON(bar_PlotsOpened_type, VALUE='Create new window', /CHECKED_MENU,	UVALUE='PLOTSOPENED_TYPE', UNAME='PLOTSOPENED_TYPE_1')
		bttn_PlotsOpened_type2 = WIDGET_BUTTON(bar_PlotsOpened_type, VALUE='Replace last window',  /CHECKED_MENU, UVALUE='PLOTSOPENED_TYPE', UNAME='PLOTSOPENED_TYPE_2')

Helpbar = WIDGET_BUTTON(barBase, VALUE='Help', UVALUE='BABOUT', /MENU)
	;bttn_usermanual = WIDGET_BUTTON(Helpbar, VALUE='User manual', UVALUE='BOPEN_USERMANUAL', SENSITIVE=0)
	bttn_about = WIDGET_BUTTON(Helpbar, VALUE='About DCE@urLAB...', UVALUE='BABOUT')

bttn_drawROI    = WIDGET_BUTTON(base_b1b1a, VALUE= 'ROI selection',  UVALUE='ROI_SELECTION')
bttn_NewROI     = WIDGET_BUTTON(base_b1b1b, VALUE= 'New ROI',  UVALUE='NEW_ROI', UNAME='NEW_ROI')

strarr_typeROIs = ['Box', 'Full', 'free']
dropl_typeROI   = WIDGET_DROPLIST(base_b1b1c, TITLE='ROI type', UVALUE='ROI_TYPE', UNAME='ROI_TYPE', $
	VALUE=strarr_typeROIs)

label_ROI_px    = WIDGET_LABEL(base_b1b2a, VALUE='Xi:')
text_ROI_px     = WIDGET_TEXT(base_b1b2a, VALUE= '0',  UVALUE='TEXT_ROI', /EDITABLE, XSIZE=3, UNAME= 'POINT_ROI_X')
label_ROI_py    = WIDGET_LABEL(base_b1b2a, VALUE='Yi:')
text_ROI_py     = WIDGET_TEXT(base_b1b2a, VALUE= '0',  UVALUE='TEXT_ROI', /EDITABLE, XSIZE=3, UNAME ='POINT_ROI_Y')
base_sep        = WIDGET_BASE(base_b1b2a, XSIZE=5)

label_ROI_size = WIDGET_LABEL(base_b1b2a, VALUE='ROI size ')
text_ROI_size_x = WIDGET_TEXT(base_b1b2a, VALUE=STRTRIM(stp.size_roi[0],2), UVALUE='TEXT_ROI',$
	/EDITABLE, XSIZE=3, UNAME= 'ROI_SIZE_X', ALL_EVENTS=0)
text_ROI_size_y = WIDGET_TEXT(base_b1b2a, VALUE=STRTRIM(stp.size_roi[1],2), UVALUE='TEXT_ROI',$
	/EDITABLE, XSIZE=3, UNAME= 'ROI_SIZE_Y', ALL_EVENTS=0)

dropl_ROIresolution_array = ['1','2','3','4','5','8','10']
dropl_ROIresolution = WIDGET_DROPLIST(base_b1b2b, TITLE='Resolution', UVALUE='TEXT_ROI',UNAME='ROI_RESOLUTION', VALUE=dropl_ROIresolution_array)

bttn_AnalyzeROI  = WIDGET_BUTTON(base_b1c1,   VALUE= 'Analyze ROI (voxel-wise)',  UVALUE='ANALYZE_ROI', UNAME='COMPUTE_ROI')

st_models_temp   = DCEMRI_ConfigModelParameters(0, ERROR_ID=error_id, ARR_MODELS=arr_models, STRARR_MODELS=strarr_models)
strarr_plasmacurve = ['Not selected', 'Biexponential', 'From File', 'From ROI']
strarr_rrcurve = ['Not selected',     'Loaded from disk', 'From ROI']

dropls_model     = WIDGET_DROPLIST(base_b1c2, TITLE='Pharm. Model ', UNAME='DROPLIST_MODEL',UVALUE='CHOOSE_MODEL', VALUE=strarr_models)

dropls_AIFcurve   = WIDGET_DROPLIST(base_b1c3, TITLE='AIF ', UNAME='DROPLIST_PLASMACURVE' ,UVALUE='CHOOSE_PLASMACURVE', VALUE=strarr_plasmacurve)
dropls_RRCurve   = WIDGET_DROPLIST(base_b1c4, TITLE='RR curve ', UNAME='DROPLIST_RRCURVE' ,UVALUE='CHOOSE_RRCURVE', VALUE=strarr_rrcurve)

WIDGET_CONTROL, dropls_model, SET_DROPLIST_SELECT=0
WIDGET_CONTROL, dropls_AIFcurve, SET_DROPLIST_SELECT=0
WIDGET_CONTROL, dropls_RRcurve, SET_DROPLIST_SELECT=0

bttn_AIFfromROI = WIDGET_BUTTON(base_b1c3, VALUE='ROI->AIF', UNAME='AIF_FROM_ROI', UVALUE='SELECT_ROI')
bttn_RRfromROI  = WIDGET_BUTTON(base_b1c4, VALUE='ROI->RR', UNAME='RR_FROM_ROI', UVALUE='SELECT_ROI')

base_rell  = WIDGET_BASE(base_b1c3, XSIZE=5)
base_rell  = WIDGET_BASE(base_b1c4, XSIZE=5)

bttn_ViewAIF = WIDGET_BUTTON(base_b1c3, VALUE='Plot', UNAME='PLOT_AIF', UVALUE='PLOT_AIF', SENSITIVE=0)
bttn_ViewRR  = WIDGET_BUTTON(base_b1c4, VALUE='Plot',  UNAME='PLOT_RR', UVALUE='PLOT_RR', SENSITIVE=0)

st_widgets_zoom    = fwBaseStruct_Zoom(PARENT_BASE=base_b1a,MAIN_BASE=base_b1a, FRAME=1)

st_widgets_param   = fwBaseStruct_Parameters(PARENT_BASE=base_b1d,	ST_PARAM=st_param, ST_AIF=st_aif)

st_widgets_frames  = fwBaseStruct_Frames(PARENT_BASE=base_b2a)

st_drawview  = fwBaseStruct_Window(PARENT_BASE=base_b2b, TITLE='Dynamic image', $
	ST_OPTIONS=st_draw_options, SIZE_IMAGE=st_inf.size_im, PALETTE_ROI=idx_paletteROI_ini, PALETTE_IMAGE=0, DRAW_ID=draw1, UVALUE='DRAW', UNAME='WINDOW1', /ROI)

strarr_text1 = STRARR(1) + '.'
text_info1   = WIDGET_TEXT(base_b2c, VALUE=strarr_text1, /ALIGN_LEFT, XSIZE=30, YSIZE=1);, FONT=font1)

st_widgets_scales = fwBaseStruct_scales(PARENT_BASE=base_b3a)

st_RCEview   = fwBaseStruct_Window(PARENT_BASE=base_b3b, TITLE='Relative contrast enhancement (RCE)', $
	ST_OPTIONS=st_draw_options, SIZE_IMAGE=st_inf.size_im, PALETTE_IMAGE=idx_paletteRCE_ini, DRAW_ID=draw2, UVALUE='DRAW', UNAME='WINDOW2')

text_info2  = WIDGET_TEXT(base_b3c, VALUE=strarr_text1, /ALIGN_LEFT, XSIZE=30, YSIZE=1);, FONT=font1)

;----------------------------------------------------------------------
st_widgets = { $
	wAuxMainBases : ULONARR(10),    $ ; reserved for children main bases
    wbase_main  : base_main, $
    wbase_smain : base_smain,$
    wbase_newROI    : base_b1b1a,   $
    wbase_typeROI   : base_b1b1c,   $
    wbase_paramROI  : base_b1b2a,   $
    wbase_resolROI  : base_b1b2b,   $
    wbase_analysis    : base_b1c1,  $
    wbase_parameters  : base_b1d,   $
    wbases_images     : [base_b2, base_b3], $
    ;---------------------------------------
    wdropls_model : dropls_model,   $
  	wdropls_AIFcurve: dropls_AIFcurve,$
	wdropls_RRCurve : dropls_RRCurve,$
	;---------------------------------------
    wfilebar : filebar,$
    wquitbar : quitbar,$
    wExportBar : exportbar,$
    wbttn_exportROIkinetics : bttn_exportROIkinetics,$
    wbttn_exportROI : bttn_exportROI,$
    wbttn_importROI : bttn_importROI,$
    woptionsbar    : optionsbar,$
    wviewbar       : viewbar, $
    wviewbar_rce   : viewbar_rce, $
    wviewbar_t1map : viewbar_t1map,$
    wviewbar_data1  : viewbar_data1,$
    wviewbar_data2  : viewbar_data2,$

	wbttns_viewbarROIs   : bttns_viewbarROIs,$
    wdropl_typeROI   : dropl_typeROI,$
    ;-----------------------------------------------
    wbttn_interp : st_widgets_zoom.bttn_interp,$
    ;-----------------------------------------------
    wslide_n     : st_widgets_frames.slide_n,$
    wslide_z     : st_widgets_frames.slide_z,$
    wlabel_n     : st_widgets_frames.label_n,$
    wlabel_z     : st_widgets_frames.label_z,$
    ;-----------------------------------------------
    wtext_minval    : st_widgets_scales.text_minval,$
    wtext_maxval    : st_widgets_scales.text_maxval,$
    wslide_minval    : st_widgets_scales.slide_minval,$
    wslide_maxval    : st_widgets_scales.slide_maxval,$
    wlabel_threshold : st_widgets_scales.label_threshold,$
    wslide_threshold : st_widgets_scales.slide_threshold,$
    ;-----------------------------------------------
    wtext_info1 : text_info1,$
    wtext_info2 : text_info2,$
    wtext_roi_px : text_roi_px,$
    wtext_roi_py : text_roi_py,$
    wtext_roi_size : [text_roi_size_x,text_roi_size_y],$
    wdropl_ROIresolution : dropl_ROIresolution,$
    wdropl_ROIresolution_array : dropl_ROIresolution_array,$
    wdraw1 : draw1,$
    wdraw2 : draw2,$
    wbttn_newROI  : bttn_newROI,$
    wbttn_drawROI : bttn_drawROI,       $
    wbttn_AnalyzeROI    : bttn_AnalyzeROI, $
    wopenbar_t1         : openbar_t1,$
    wopenbar_t1map      : openbar_t1map,$
    wopenbar_t1free     : openbar_t1free, $
    wopenbar_t1mapfree  : openbar_t1mapfree, $
    wopenbar_AIF        : openbar_AIF,$
    wopenbar_RR         : openbar_RR, $
    wbttn_viewAIF       : bttn_viewAIF, $
    wbttn_viewRR        : bttn_viewRR,  $
    wbarmenu_limitframes : barmenu_limitframes,$
	wbttn_RCE_type1 : bttn_RCE_type1,$
	wbttn_RCE_type2 : bttn_RCE_type2,$
	wbttn_PlotsOpened : [bttn_PlotsOpened_type1, bttn_PlotsOpened_type2], $
    last : 0}
;----------------------------------------------------------------------

n_images_right = 3; (RCE, IAUC y TTM)
st_roi             = fStruct_Create(/ROISTRUCT)
st_roi.Np_lr       = stp.np_lr
st_roi.opt_typeroi = stp.type_roi
st_results         = fStruct_Create(/RESULTSTRUCT)

St={ $
	data   : PTR_NEW(/ALLOCATE_HEAP),$             ; dynamic image DCE-MRI
    mask   : PTR_NEW(/ALLOCATE_HEAP),$             ; auxiliar mak
    t1map  : PTR_NEW(/ALLOCATE_HEAP),$             ; T10 map
    parametric_AIF  : PTR_NEW(/ALLOCATE_HEAP),$    ; parametric values of AIF (or Cp) are stored in a pointer, because the model can vary
    modelinfo       : PTR_NEW(/ALLOCATE_HEAP), $
	inf       : st_inf,          $
	opt       : st_draw_options, $
	par       : st_param,        $
	ptrarr_r  : PTRARR(20, /ALLOCATE_HEAP),   $ ; Whe load until 20 results (NOT IMPLEMENTED YET)
	results   : st_results,      $
	roi       : st_roi,     $
    wd  : st_widgets,       $      ; Widget identifiers
    wdp : st_widgets_param, $      ; Widget identifiers (parameters)
    vi1 : st_drawview,  $          ; Graphic objects od display 1
    vi2 : st_RCEview,   $          ; Graphic objects od display 2 (RCE)
    objplots : PTR_NEW(/ALLOCATE_HEAP), $ ; Pointer to configuration of graphic objects (read from .sav file)

	flag_stddev : 0l, $
	flags_applyparams    : [1,1,0L,0L,0L], $
	flag_drawROI        : 0l,$
	flag_analysis       : 0l,$
	flag_analysis_done  : 0l,$
	flag_interp         : stp.flag_interpolation,$
	flag_send_Event     : 0l,$
	flag_noresetroi     : 0l,$
	flag_hideroilimits : 0l, $
	flag_freeROIclose  : 0L, $
	flag_freeROIinit   : 0L, $
	flag_freeROIcontinuous : 0l, $
	flag_aif_from_roi      : 0l, $
	flag_aif_from_file     : 0l, $
	flag_RR_from_roi       : 0l, $
	flag_RR_from_file      : 0l, $

    opt_drawROIoutliers  : 2l, $ ; ROI outlayers (outside of limits of SCALE) are: 0) drawn 1) transparent 2) black
    opt_typeRCE        : 1l, $ ; tipo de calculo RCE (1, relativo a valor medio antes de contraste:  2 - Absoluto max-min)
    opt_current_param  : 0l, $ ; actual visualized parameters
    opt_current_model  : 0l, $ ; actual model
    opt_AIFcurve       : 0l, $ ; Tipo de curva de CA en plasma usado (1- Biexponencial analítica, 2 - de Archivo, 3 - de ROI)
    opt_RRcurve        : 0l, $ ; Tipo de curva de RR usada    (1 - de Archivo, 2 - de ROI)
    opt_view_im_right  : 0l, $ ; opción de visualización a la derecha (RCE=0, AIUC=1, TTM=2)
    opt_plotsOpened    : 2l, $ ; en la opcion por defecto, las ventanas machacan a las anteriores (2) en la otra, se van poniendo encima
    initial_zoom   : initial_zoom,   $  ; Zoom sobre [256x256]
    path_data      : stp.path_data,  $
    path_t1map     : stp.path_t1map, $
    path_result    : stp.path_result,$
    path_aif       : stp.path_aif,   $
    opt_nodestroyatend : 0B, $
    changing_palette : 0, $

    wtam_z: 0l $
}
pSt = PTR_NEW(St, /NO_COPY)

IF SIZE(st_aif, /TNAME) EQ 'STRUCT' THEN BEGIN
	*(*Pst).parametric_AIF = st_aif
ENDIF ELSE BEGIN
	*(*Pst).parametric_AIF = -1 ; at least one value
ENDELSE

size_b3 = (WIDGET_INFO(base_b, /GEOMETRY)).scr_xsize
size_b1 = (WIDGET_INFO(base_b1, /GEOMETRY)).scr_xsize

size_temp = (WIDGET_INFO(base_b3a, /GEOMETRY)).scr_ysize
WIDGET_CONTROL, base_b2a, SCR_YSIZE=size_temp

;--------------------------------------------------------------
OPT_view_RCE = 1
IF opt_view_RCE EQ 0 THEN BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbases_images[1], MAP=0, XSIZE=1, YSIZE=1
ENDIF ELSE BEGIN
	WIDGET_CONTROL, (*Pst).wd.wbases_images[1], MAP=1
ENDELSE
WIDGET_CONTROL, (*Pst).wd.wviewbar_rce, SET_BUTTON=opt_view_RCE
;--------------------------------------------------------------

size_button_y = (WIDGET_INFO(bttn_drawROI, /GEOMETRY)).scr_ysize

ok = SetWidget_maxsizes([base_b1a, base_b1b, base_b1c, base_b1d], /XSIZE)
ok = SetWidget_maxsizes([dropls_model,dropls_AIFcurve,dropls_RRCurve], /XSIZE)

pos_resolution = (WHERE((*Pst).wd.wdropl_ROIresolution_Array EQ stp.np_lr,ct))[0]
IF ct NE 1 THEN BEGIN
	ok = DIALOG_MESSAGE('ROI Resolution = ' + STRTRIM(stp.np_lr,2)+ ' not supported, change to 1', /INFO)
	(*Pst).roi.np_lr = 1L
ENDIF

ok = fwdcemri_activate_bases(Pst, /INITIAL_STATE)

IF (*Pst).opt_typeRCE EQ 1 THEN WIDGET_CONTROL, (*Pst).wd.wbttn_RCE_type1, SET_BUTTON=1
IF (*Pst).opt_typeRCE EQ 2 THEN WIDGET_CONTROL, (*Pst).wd.wbttn_RCE_type2, SET_BUTTON=1
IF (*Pst).opt_plotsOpened EQ 1 THEN WIDGET_CONTROL, (*Pst).wd.wbttn_plotsOpened[0], SET_BUTTON=1
IF (*Pst).opt_plotsOpened EQ 2 THEN WIDGET_CONTROL, (*Pst).wd.wbttn_plotsOpened[1], SET_BUTTON=1

ok = GetWidget_PositionInScreen((*Pst).wd.wbase_main, CENTER=center)

LOADCT, 0

;-------------------------------------------------------------------
file_configplots = '.\Config\Configplots.pref'
objplots = Restore_GraphicObjectsData(FILE=file_configplots)
IF SIZE(objplots,/TNAME) EQ 'STRUCT' THEN BEGIN
	*(*Pst).objplots = objplots
ENDIF ELSE *(*Pst).objplots = -1
;-------------------------------------------------------------------

WIDGET_CONTROL, (*Pst).wd.wbase_main, /REALIZE, /MAP, SET_UVALUE=pSt

;-----------------------------------------------------------------------
XMANAGER, 'Interface_DCEMRI', (*Pst).wd.wbase_main, NO_BLOCK=1, $
    EVENT_HANDLER='Interface_DCEMRI_event', CLEANUP='Interface_DCEMRI_Destroy'
;-----------------------------------------------------------------------

HELP, /MEMORY

END ; Interface_DCEMRI

;*************************************************************************************************
;*************************************************************************************************
;<-->
;<++>

FUNCTION fwidget_Plot_Model_v2, Pst,  WIN_NUMBER=win_number, POINT=point, $
	P_IDX=p_idx, N_HD=n_hd, ARR_DATA=arr_data, ST_INFO=st_info, ARR_MODEL=arr_model, ARR_SIGNAL=arr_signal

; N_HD -> Number of interpolation point between samples for high definition plot

IF N_ELEMENTS(n_hd) EQ 0 THEN n_hd=1L

;(*(*Pst).modelinfo).using_Ct  ; Plots Ct or Signal intensity..

IF N_ELEMENTS(p_idx) NE 0  THEN opt_internal = 1 ELSE opt_internal = 0

IF N_ELEMENTS(point) EQ 2 THEN BEGIN
	 str_point_lr = ' at ROI point: [' + STRTRIM(point[0],2) + ', ' + STRTRIM(point[1],2) + ']
ENDIF ELSE str_point_lr = ''

CASE (*(*Pst).modelinfo).using_Ct OF
0 : BEGIN  ; the model plots the signal (RCE relative to the initial value. 100%)
	IF opt_internal EQ 0 THEN BEGIN
		arr_plot  = REFORM(arr_data)
	ENDIF ELSE BEGIN
		;arr_Signal = REFORM((*(*Pst).results.data_Signal_lr)[p_idx,*])
		;arr_RCE = RCE_signalRel_from_Signal(arr_Signal, FRAME_INJECTION=(*Pst).par.nframe_injection)
		arr_RCE = REFORM((*(*Pst).results.data_signal_RCE)[p_idx,*])
		arr_plot= arr_RCE*100
		undefine, arr_RCE, arr_signal
	ENDELSE
	ytitle = 'S(t)'
	title  = 'MRI signal (%) relative to initial before CA injection)' + str_point_lr
END
1 : BEGIN ; the model plots the Ct concentration
	IF opt_internal EQ 0 THEN BEGIN
		arr_plot  = REFORM(arr_data)
	ENDIF ELSE BEGIN
		arr_ct   = REFORM((*(*Pst).results.data_CT)[p_idx,*])
		arr_plot = TEMPORARY(arr_ct)
	ENDELSE
	ytitle = 'C(t)'
	title  = 'Tissue concentration (mM or mmol/litre)' + str_point_lr
END
ENDCASE

arr_t_min        = *(*Pst).results.arr_time
parameters_ids  = (*(*Pst).modelinfo).ids
parameters_name = (*(*Pst).modelinfo).idsinf[2,*]

n_parameters = N_ELEMENTS(parameters_ids)

IF N_ELEMENTS(parameters) EQ 0 THEN BEGIN
	parameters = FLTARR(n_parameters)
    standard_dev = FLTARR(n_parameters)
	FOR i=0l, n_parameters-1 DO BEGIN
		parameters[i]     =  (*(*Pst).results.dataptr[i,0])[p_idx]
		IF (*(*Pst).modelinfo).sd_valid[i] THEN $
			standard_dev[i]   =  (*(*Pst).results.dataptr[i,1])[p_idx]
	ENDFOR
ENDIF

IF (n_parameters EQ N_ELEMENTS(parameters)) EQ 0 THEN RETURN,-1

xtitle =     'Time (min)'

st_info = {parameters:parameters, $
		   standard_dev:standard_dev,$
		   ytitle:ytitle, $
		   xtitle:xtitle, $
		   title :title,  $
		   n_hd : n_hd}


;-------------------------------------------------------------------
arr_t_estimated = INDGEN(N_ELEMENTS(arr_t_min)*n_hd-(n_hd-1))*1.0
arr_t_estimated/= MAX(arr_t_estimated)
arr_t_estimated*= (MAX(arr_t_min)-arr_t_min[0]) + arr_t_min[0]
;-------------------------------------------------------------------
time_baseline    = arr_t_min[(*Pst).par.nframe_injection]
arr_t_est_adjust = arr_t_estimated[(*Pst).par.nframe_injection*n_hd: N_ELEMENTS(arr_t_estimated)-1]
arr_t_est_adjust-= arr_t_estimated[(*Pst).par.nframe_injection*n_hd]
;-------------------------------------------------------------------

opt_model = 1
CASE (*(*Pst).modelinfo).typemodel OF

3  : BEGIN  ; Tofts o Tofts extendido

	CASE (*Pst).results.tofts_type OF
	1 : BEGIN ; convolution
		arr_plot_estimated = Model_Ct_Tofts_convol(ARR_TIME=arr_t_estimated, ARR_CP=*(*Pst).results.arr_cp, PARAMS=parameters, NORMALIZE=normalize, $
				METHOD=1, IDS_INPUT=parameters_ids) ; method 1 is convolution
		arr_ct_estimated = arr_plot_estimated
	END
	2 : BEGIN ; analitic
		dose    = (*Pst).par.dose
		cpm    = *(*Pst).results.cp_model
		model_aif = [cpm.a1, cpm.a2, cpm.m1, cpm.m2]	; only biexponential model used at this time

		arr_ct_estimated   = Model_Ct_Tofts(DOSE=dose, ARR_TIME=arr_t_est_adjust, PARAMS=parameters, MODEL_AIF=model_aif, IDS_INPUT=parameters_ids)
		arr_plot_estimated = [REPLICATE(0, (*Pst).par.nframe_injection*n_hd), arr_ct_estimated]
	END
	ENDCASE
END
2  : BEGIN  ; hoffmann
	arr_Rce_estimated  = Model_RCE_Hoffmann(ARR_TIME=arr_t_est_adjust, PARAMS=parameters, IDS_INPUT=parameters_ids)
	arr_Rce_estimated*=100
	IF (*Pst).par.nframe_injection GT 0 THEN $
		arr_plot_estimated = [REPLICATE(100, (*Pst).par.nframe_injection*n_hd), arr_Rce_estimated] $
	ELSE arr_plot_estimated = arr_RCE_estimated

END
4  : BEGIN  ; Larsson

	model_aif = [(*Pst).par.a1, (*Pst).par.a2,(*Pst).par.m1, (*Pst).par.m2]
	arr_Rce_estimated = Model_RCE_Larsson(ARR_TIME=arr_t_est_adjust, PARAMS=parameters, MODEL_AIF=model_aif, IDS_INPUT=parameters_ids)
	arr_Rce_estimated*=100
	IF (*Pst).par.nframe_injection GT 0 THEN $
		arr_plot_estimated = [REPLICATE(100, (*Pst).par.nframe_injection*n_hd), arr_Rce_estimated] $
	ELSE arr_plot_estimated = arr_RCE_estimated

END
5  : BEGIN ; RR moswl
	arr_ct_RR          = *(*Pst).results.arr_ct_rr
	arr_plot_estimated = Model_Ct_Yankeelov_convol(ARR_TIME=arr_t_estimated, ARR_CT_RR=arr_ct_RR, PARAMS=parameters, IDS_INPUT=parameters_ids,$
							KTRANS_RR=(*Pst).par.ktrans_rr, KEP_RR=(*Pst).par.kep_rr, NORMALIZE=normalize, METHOD=1)
	;arr_plot_estimated = [REPLICATE(0, (*Pst).par.nframe_injection*n_hd), arr_ct_estimated]
	arr_ct_estimated = arr_plot_estimated
END
1  : BEGIN  ; no model
	opt_model = 0L
	arr_plot_estimated = arr_plot
END
ENDCASE
IF ((N_ELEMENTS(arr_RCE_estimated) LE 1) AND (N_ELEMENTS(arr_CT_estimated) LE 1) AND (opt_model EQ 1)) THEN BEGIN
	ok = DIALOG_MESSAGE('Critical error', /INFO) & RETURN,-1
ENDIF
undefine, arr_ct_estimated, arr_RCE_estimated

arr_signal = [[arr_plot],[arr_t_min]]

IF opt_model EQ 1 THEN BEGIN
	arr_model  = [[arr_plot_estimated],[arr_t_estimated]]
ENDIF ELSE BEGIN
	undefine, arr_model, arr_plot_estimated ,arr_t_estimated
ENDELSE

RETURN ,1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemrirs_AdjustScales, Pst

ni = (SIZE((*Pst).results.dataptr, /DIMENSIONS))[0]
nj = (SIZE((*Pst).results.dataptr, /DIMENSIONS))[1]

n_points_lr = N_ELEMENTS(*(*Pst).roi.idx_points_lr)

FOR j=0l, nj-1 DO BEGIN
	FOR i=0l, ni-1 DO BEGIN

	IF 	SIZE(*(*Pst).results.dataptr[i,j], /TYPE) NE 0 THEN BEGIN ; no está indefinido
		min_val = MIN((*(*Pst).results.dataptr[i,j])[0:n_points_lr-1], MAX=max_val)
		(*Pst).results.scales[*,i,j] = [min_val, max_val]
	ENDIF

	ENDFOR
ENDFOR

RETURN, 1

END

;**************************************************************************************************
;**************************************************************************************************

PRO fwdcemrirs_Change_Palette1b, DATA=data

    TVLCT, R,G,B, /GET

    data.vi1.opalette    -> setproperty, RED=r,GREEN=G, BLUE=b
    ;data.vi1.ohcolorbar -> setproperty, PALETTE  = data.vi1.opalette

    WIDGET_CONTROL, data.wd.wDraw1, GET_VALUE = window_1
    window_1 -> Draw, data.vi1.oView

END

;**************************************************************************************************
;**************************************************************************************************

PRO fwdcemrirs_Change_Palette2b, DATA=data ; paleta de color de la roi

    TVLCT, R,G,B, /GET
    data.vi1.opaletteROI   -> setproperty, RED=r,GREEN=G, BLUE=b
    data.vi1.ohcolorbar    -> setproperty, PALETTE  = data.vi1.opaletteROI

    WIDGET_CONTROL, data.wd.wDraw1, GET_VALUE = window_1
    window_1 -> Draw, data.vi1.oView

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemrirs_Plot_Model, Pst,  WIN_NUMBER=win_number, POINT=point, $
	P_IDX=p_idx, XPOS=xpos, YPOS= ypos, N_HD=n_hd, $
	ARR_DATA=arr_data, PARAMETERS = parameters, STANDARD_DEV=standard_dev

; N_HD -> Number of interpolation point between samples for high definition plot

IF N_ELEMENTS(n_hd) EQ 0 THEN n_hd=1L

;(*(*Pst).modelinfo).using_Ct  ; Plots Ct or Signal intensity..

IF N_ELEMENTS(p_idx) NE 0  THEN opt_internal = 1 ELSE opt_internal = 0

IF N_ELEMENTS(point) EQ 2 THEN BEGIN
	 str_point_lr = ' at ROI point: [' + STRTRIM(point[0],2) + ', ' + STRTRIM(point[1],2) + ']
ENDIF ELSE str_point_lr = ''

CASE (*(*Pst).modelinfo).using_Ct OF
0 : BEGIN  ; the model plots the signal (RCE relative to the initial value. 100%)
	IF opt_internal EQ 0 THEN BEGIN
		arr_plot  = REFORM(arr_data)
	ENDIF ELSE BEGIN
		;arr_Signal = REFORM((*(*Pst).results.data_Signal_lr)[p_idx,*])
		;arr_RCE = RCE_signalRel_from_Signal(arr_Signal, FRAME_INJECTION=(*Pst).par.nframe_injection)
		arr_RCE = REFORM((*(*Pst).results.data_signal_RCE)[p_idx,*])
		arr_plot= arr_RCE*100
		undefine, arr_RCE, arr_signal
	ENDELSE
	ytitle = 'S_t'
	title  = 'MRI signal (%) relative to initial before CA injection)' + str_point_lr
END
1 : BEGIN ; the model plots the Ct concentration
	IF opt_internal EQ 0 THEN BEGIN
		arr_plot  = REFORM(arr_data)
	ENDIF ELSE BEGIN
		arr_ct   = REFORM((*(*Pst).results.data_CT)[p_idx,*])
		arr_plot = TEMPORARY(arr_ct)
	ENDELSE
	ytitle = 'C_t'
	title  = 'Tissue concentration (mM or mmol/litre)' + str_point_lr
END
ENDCASE

arr_t_min = *(*Pst).results.arr_time
parameters_ids  = (*(*Pst).modelinfo).ids
parameters_name = (*(*Pst).modelinfo).idsinf[2,*]

n_parameters = N_ELEMENTS(parameters_ids)

IF N_ELEMENTS(parameters) EQ 0 THEN BEGIN
	parameters = FLTARR(n_parameters)
    standard_dev = FLTARR(n_parameters)
	FOR i=0l, n_parameters-1 DO BEGIN
		parameters[i]     =  (*(*Pst).results.dataptr[i,0])[p_idx]
		IF (*(*Pst).modelinfo).sd_valid[i] THEN $
			standard_dev[i]   =  (*(*Pst).results.dataptr[i,1])[p_idx]
	ENDFOR
ENDIF

IF (n_parameters EQ N_ELEMENTS(parameters)) EQ 0 THEN RETURN,-1

;-------------------------------------------------------------------
arr_t_estimated = INDGEN(N_ELEMENTS(arr_t_min)*n_hd-(n_hd-1))*1.0
arr_t_estimated/= MAX(arr_t_estimated)
arr_t_estimated*= (MAX(arr_t_min)-arr_t_min[0]) + arr_t_min[0]
;-------------------------------------------------------------------
time_baseline    = arr_t_min[(*Pst).par.nframe_injection]
arr_t_est_adjust = arr_t_estimated[(*Pst).par.nframe_injection*n_hd: N_ELEMENTS(arr_t_estimated)-1]
arr_t_est_adjust-= arr_t_estimated[(*Pst).par.nframe_injection*n_hd]
;-------------------------------------------------------------------

opt_model = 1
CASE (*(*Pst).modelinfo).typemodel OF

3  : BEGIN  ; Tofts o Tofts extendido

	CASE (*Pst).results.tofts_type OF
	1 : BEGIN ; convolution
		arr_plot_estimated = Model_Ct_Tofts_convol(ARR_TIME=arr_t_estimated, ARR_CP=*(*Pst).results.arr_cp, PARAMS=parameters, NORMALIZE=normalize, $
				METHOD=1, IDS_INPUT=parameters_ids) ; method 1 is convolution
		arr_ct_estimated = arr_plot_estimated
	END
	2 : BEGIN ; analitic
		dose    = (*Pst).par.dose
		cpm    = *(*Pst).results.cp_model
		model_aif = [cpm.a1, cpm.a2, cpm.m1, cpm.m2]	; only biexponential model used at this time

		arr_ct_estimated   = Model_Ct_Tofts(DOSE=dose, ARR_TIME=arr_t_est_adjust, PARAMS=parameters, MODEL_AIF=model_aif, IDS_INPUT=parameters_ids)
		IF (*Pst).par.nframe_injection GT 0 THEN $
			arr_plot_estimated = [REPLICATE(0, (*Pst).par.nframe_injection*n_hd), arr_ct_estimated] $
	 	ELSE arr_plot_estimated = arr_ct_estimated
	END
	ENDCASE
END
2  : BEGIN  ; hoffmann
	arr_Rce_estimated  = Model_RCE_Hoffmann(ARR_TIME=arr_t_est_adjust, PARAMS=parameters, IDS_INPUT=parameters_ids)
	arr_Rce_estimated*=100
	IF (*Pst).par.nframe_injection GT 0 THEN $
		arr_plot_estimated = [REPLICATE(100, (*Pst).par.nframe_injection*n_hd), arr_Rce_estimated] $
	ELSE arr_plot_estimated = arr_rce_estimated
END
4  : BEGIN  ; Larsson

	model_aif = [(*Pst).par.a1, (*Pst).par.a2,(*Pst).par.m1, (*Pst).par.m2]
	arr_Rce_estimated = Model_RCE_Larsson(ARR_TIME=arr_t_est_adjust, PARAMS=parameters, MODEL_AIF=model_aif, IDS_INPUT=parameters_ids)
	arr_Rce_estimated*=100
	IF (*Pst).par.nframe_injection GT 0 THEN $
	arr_plot_estimated = [REPLICATE(100, (*Pst).par.nframe_injection*n_hd), arr_Rce_estimated] $
	ELSE arr_plot_estimated = arr_rce_estimated
END
5  : BEGIN ; RR moswl
	arr_ct_RR          = *(*Pst).results.arr_ct_rr
	arr_plot_estimated = Model_Ct_Yankeelov_convol(ARR_TIME=arr_t_estimated, ARR_CT_RR=arr_ct_RR, PARAMS=parameters, IDS_INPUT=parameters_ids,$
							KTRANS_RR=(*Pst).par.ktrans_rr, KEP_RR=(*Pst).par.kep_rr, NORMALIZE=normalize, METHOD=1)
	;arr_plot_estimated = [REPLICATE(0, (*Pst).par.nframe_injection*n_hd), arr_ct_estimated]
	arr_ct_estimated = arr_plot_estimated
END
1  : BEGIN
	opt_model = 0L
	arr_plot_estimated = arr_plot
END
ENDCASE
IF ((N_ELEMENTS(arr_RCE_estimated) LE 1) AND (N_ELEMENTS(arr_CT_estimated) LE 1) AND (opt_model EQ 1)) THEN BEGIN
	ok = DIALOG_MESSAGE('Critical error', /INFO) & RETURN,-1
ENDIF
undefine, arr_ct_estimated, arr_RCE_estimated


xsize = 900
ysize = 500

IF N_ELEMENTS(xpos) EQ 0 THEN xpos=20
IF N_ELEMENTS(ypos) EQ 0 THEN ypos=20

colors = LONARR(2)
colors[0] = Get_color([200,0,0], /HEX)    ; rojo
colors[1] = Get_color([56,56,200], /HEX)  ; azul bonito

IF (*(*Pst).modelinfo).using_Ct THEN BEGIN
	minplot_abs_max =  0.0 ; value in concentration

ENDIF ELSE BEGIN
	minplot_abs_max = 95.0 ; VALUE in %
ENDELSE

diff_x = (arr_t_min[1]-arr_t_min[0])
xrange = [MIN([arr_t_min,arr_t_estimated])-diff_x/2.0, MAX([arr_t_min,arr_t_estimated])+diff_x/2.0]
yrange = [MIN([arr_plot,arr_plot_estimated]) < minplot_abs_max, MAX([arr_plot,arr_plot_estimated])*1.1]

xtitle = 'Time (min)'
str_set_font = "Arial*16"

type_data_estimated = 1 ; continuous line
IF N_ELEMENTS(arr_plot) GT 100 THEN type_data_ct=1 ELSE type_data_ct=2
;---------------------------------------------------------------------------------
ok = PlotProfile(DATA_X=arr_t_min, DATA_Y=arr_plot, WIN_NUMBER=win_number, $
		TITLE=title, XTITLE=xtitle, YTITLE=ytitle, $
		YRANGE=yrange,SET_FONT=str_set_font,$
		STR_INFO=str_info, XPOS=xpos, YPOS=ypos, XSIZE=xsize,YSIZE=ysize,$
		TYPE_DATA=type_data_ct, COLORS=colors[0], THICKS=thicks, COLOR_BACKGROUND='FFFFFF'x, $
		POSITION=[0.06, 0.1, 0.70, 0.95])

IF opt_model EQ 1 THEN BEGIN
	;---------------------------------------------------------------------------------
	ok = PlotProfile(DATA_X=arr_t_estimated, DATA_Y=arr_plot_estimated, OPLOTT=1, $
			SET_FONT=str_set_font, STR_INFO=str_info,$
			TYPE_DATA=type_data_estimated, COLORS=colors[1], THICKS=thicks,COLOR_BACKGROUND='FFFFFF'x)
	;---------------------------------------------------------------------------------
	;WSET, 1

	str_model = (*(*Pst).modelinfo).str_model

	nn_pos = (*(*Pst).modelinfo).nparam

	arr_text = STRARR(2,nn_pos+2)
	arr_pos  = INDGEN(nn_pos)
	arr_ip   = INDGEN(nn_pos)+1
	nn_pos = (*(*Pst).modelinfo).nparam
	arr_text[0,0] = str_model
	opt_error = 1
	FOR i=0l, nn_pos-1 DO BEGIN
		arr_text[0,arr_ip[i]] =  (*(*Pst).modelinfo).idsinf[2,arr_pos[i]] + ' : '
		arr_text[1,arr_ip[i]] = STRTRIM(STRING(parameters[arr_pos[i]],  FORMAT='(F15.4)'),2)
		opt_error = (*(*Pst).modelinfo).sd_valid[arr_pos[i]]
		IF opt_error THEN BEGIN
			IF N_ELEMENTS(standard_dev[arr_pos[i]]) GT 0 THEN BEGIN
				arr_text[1,arr_ip[i]]+= ' ± '+ STRTRIM(STRING(standard_dev[arr_pos[i]],  FORMAT='(F10.4)'),2)
			ENDIF
		ENDIF
	ENDFOR

	ok = Plot_xyouts_matrix(arr_text, XPOS=0.72, YPOS=0.90, COLOR='00'x, FONT=0, $
		SET_FONT=str_set_font, SEP_COL=1, BOLD=1)
ENDIF

END

;**************************************************************************************************
;**************************************************************************************************


FUNCTION fwdcemrirs_DrawingEvents_Result, ev, Pst

    WIDGET_CONTROL, (*Pst).wd.wDraw1,  GET_VALUE = window_draw

    possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL', 'EXPOSE' ]
    EventType = possibleEventTypes[ev.type]

	error = window_draw -> PickData((*Pst).vi1.oview, (*Pst).vi1.omodel, [ev.x,ev.y], Coord_view)
	id_text = (*Pst).wd.wtext_info1
	;------------------------------------------------------------------------------
    IF (coord_view[0] LT 0 ) OR (coord_view[1] LT 0) THEN BEGIN
    	WIDGET_CONTROL, id_text, SET_VALUE=''
    	RETURN,-1
    ENDIF
    IF (coord_view[0] GE (*Pst).inf.size_im[0]) OR (coord_view[1] GE (*Pst).inf.size_im[1]) THEN BEGIN
    	WIDGET_CONTROL, id_text, SET_VALUE=''
    	RETURN,-1
    ENDIF
    point_fix = FIX(coord_view[0:1])
	;------------------------------------------------------------------------------
	type_action = ''
	opt_fast = 0l
	;------------------------------------------------------------------------------

	CASE EventType OF
	'MOTION': BEGIN
		type_action = 'Write_position_and_value'
	END
	'DOWN' : BEGIN ; ev.press = 1, botón izquierdo, ev.press = 4, botón derecho
		type_action = 'Write_ROI_value_inMotion'
		IF (ev.press EQ 4)  THEN $
			type_action = 'Draw_params_in_ROI_box'
	END
	ELSE: RETURN,-1
	ENDCASE

	IF opt_fast EQ 0 THEN BEGIN
		WIDGET_CONTROL,(*Pst).wd.wslide_z,  GET_VALUE=n_z
    	WIDGET_CONTROL,(*Pst).wd.wslide_n,  GET_VALUE=n_n
		tam = [(*Pst).inf.size_im[0],(*Pst).inf.size_im[1],(*Pst).inf.nslices,(*Pst).inf.nframes]
	END

	CASE type_action OF
	;------------------------------------------------------------------------------
	'Write_position' : BEGIN
		point_v2 = get_point_YupConvention(point_fix, DIMENSIONS=(*Pst).inf.size_im)
		ok = fwdcemri_textWrite_pointval(id_text, point_v2)
	END
	;------------------------------------------------------------------------------
	'Write_position_and_value' : BEGIN
		val = (*(*Pst).data)[point_fix[0], point_fix[1], n_z, n_n]
		point_v2 = get_point_YupConvention(point_fix, DIMENSIONS=(*Pst).inf.size_im)
		ok = fwdcemri_textWrite_pointval(id_text, point_v2, VALUE=val, UNITS=str_units)
	END
	;------------------------------------------------------------------------------
 	'Draw_params_in_ROI_box' : BEGIN

   		n_hd = 1l
   		IF n_z NE (*Pst).roi.slicez THEN RETURN,-1

   		val = (*(*Pst).data)[point_fix[0], point_fix[1], (*Pst).roi.slicez, n_n]

		point_v2 = get_point_YupConvention(point_fix, DIMENSIONS=(*Pst).inf.size_im)
		ok = fwdcemri_textWrite_pointval(id_text, point_v2, VALUE=val, UNITS=str_units)

   		p_idx = ROImap_GetPoint(ST_ROI=(*Pst).roi, point_fix, POINT_IMF=pplr_if) ; Obtain point index in ROI
		IF p_idx EQ -1 THEN RETURN,-1

		opt_plot = 1
		CASE opt_plot OF
		1 : BEGIN
			win_number  = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)

			ok = fwdcemrirs_Plot_Model(Pst,  WIN_NUMBER=win_number, POINT=point, $
				P_IDX=p_idx, XPOS=xpos, YPOS= ypos, N_HD=n_hd)
		END
		2 : BEGIN
			ok = fwidget_Plot_Model_v2(Pst, POINT=point, P_IDX=p_idx, N_HD=n_hd, ARR_DATA=arr_data, $
				ARR_SIGNAL=arr_signal, ARR_MODEL=arr_model, ST_INFO=st_info)

			;--------------------------------------------------------------------------------------------
			Interface_Plot, ST_OPTIONS=st_options, ST_MODELINFO=*(*Pst).modelinfo, ST_INFO=st_info, $
				ARR_SIGNAL=arr_signal, ARR_MODEL=arr_model, XPOS=xpos, YPOS=ypos,$
				PARENT_BASE=(*Pst).wd.wbase_main, MAIN_BASE=(*Pst).wbase_plot, BASE_OUTPUT=base_output
			;--------------------------------------------------------------------------------------------
			(*Pst).wbase_plot =base_output
			;--------------------------------------------------------------------------------------------
		END
		ENDCASE

	END
	;------------------------------------------------------------------------------
    'draw_dynProfile_LowRes' : BEGIN
    	p_idx = ROImap_GetPoint(ST_ROI=(*Pst).roi, point_fix, POINT_IMF=pplr_if)

		win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
		IF p_idx NE -1 AND ((*Pst).roi.opt_typeROI) EQ 0 THEN BEGIN
			arr_dynMRI_lr = REFORM((*(*Pst).results.data_Signal_lr)[p_idx,*])
			ok = plot_profile_dynMRI(arr_dynMRI_lr, ARR_TIME=arr_time, XPOS=xpos, YPOS=ypos,$
				WIN_NUMBER=win_number, POINT=pplr_if, FRAME_INJECTION = (*Pst).par.nframe_injection)
		ENDIF ELSE BEGIN ; fuera de la ROI
			arr_dyn = REFORM((*(*Pst).data)[point_fix[0],point_fix[1],n_z,*])
			ok = plot_profile_dynMRI(arr_dyn, WIN_NUMBER=win_number, HR_POINT=point_fix[0:1],$
				FRAME_INJECTION = (*Pst).par.nframe_injection, XPOS=xpos, YPOS=ypos)
		ENDELSE
	END
	;------------------------------------------------------------------------------
    'draw_dynProfile_HighRes' : BEGIN
    	arr_dyn = REFORM((*(*Pst).data)[point_fix[0],point_fix[1],n_z,*])

    	opt_drawprofile = 1
    	opt_gaussian    = 0
    	;(*Pst).par.nframe_injection  = 3l
    	; ultima frame que se utiliza para calcular la concentracion sin contraste
    	IF opt_gaussian THEN $
			arr_dyn = Filter_Gaussian_1D(arr_dyn, SIGMA=3, /FLOATING)

    	win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)

    	CASE opt_drawprofile OF
    	1 : BEGIN ; direct profile
			ok = plot_profile_dynMRI(arr_dyn, ARR_TIME=*(*Pst).results.arr_time, WIN_NUMBER=win_number, HR_POINT=point_fix[0:1],$
				FRAME_INJECTION = (*Pst).par.nframe_injection, XPOS=xpos, YPOS=ypos)
		END
		2 : BEGIN ; profile of relative increment of susceptibility (DSC images) R2 = -ln(S(t)/S(0))
			;arr_dyn[0:100] = S_o
			S_o = MEAN(arr_dyn[0:(*Pst).par.nframe_injection])
			ct_rel = -ALOG(arr_dyn/S_o)
			ok = plot_profile(ct_rel, WIN_NUMBER=win_number, HR_POINT=point_fix[0:1], NO_PLOT_PSYM=0, XPOS=xpos, YPOS=ypos)
		END
		ELSE:
		ENDCASE
		opt_prueba = 0
		IF opt_prueba EQ 1 THEN BEGIN
		;----------------------------------------
			im = REFORM((*(*Pst).data)[*,*,n_z,*])
			tot_begin = TOTAL(im[*,*,0:20-1],3)/20.0
			tot_max = im[*,*,0]*0
			FOR i=0l,(*Pst).inf.nframes-1 DO BEGIN
				tot_max = tot_max > im[*,*,i]
			ENDFOR
			im_diff = tot_max/tot_begin

			viewg, im_diff GT 1.5
		ENDIF

    END
    ;------------------------------------------------------------------------------
    'Write_ROI_value_inMotion': BEGIN

		;IF (*Pst).opt_typeROI EQ 2 THEN RETURN ; solo para ROIs cuadradas

		p_idx = ROImap_GetPoint(ST_ROI=(*Pst).roi, point_fix, POINT_IMF=pplr_if)

		IF p_idx EQ -1 THEN BEGIN
			point_v2 = get_point_YupConvention(point_fix, DIMENSIONS=(*Pst).inf.size_im)
			ok = fwdcemri_textWrite_pointval(id_text, point_v2)
		ENDIF ELSE BEGIN
			IF N_ELEMENTS(pplr_if) EQ 0 THEN RETURN,-1

			val_str = ''

			i_stddev = (*Pst).flag_stddev*(*(*Pst).modelinfo).sd_valid[(*Pst).opt_current_param] ;se escribe o no la desviación estándar
			val = (*(*Pst).results.dataptr[(*Pst).opt_current_param,i_stddev])[p_idx]
			val_str = (*(*Pst).modelinfo).idsinf[0,(*Pst).opt_current_param] + STRTRIM(STRING(val, FORMAT='(F8.2)'),2) +   (*(*Pst).modelinfo).idsinf[1,(*Pst).opt_current_param]
			IF i_stddev EQ 1 THEN val_str+=' (SD)'

			WIDGET_CONTROL, id_text, SET_VALUE= ' ROI: [' + STRTRIM(pplr_if[0],2) + ', ' + $
			 	STRTRIM(pplr_if[1],2) + '] ' + val_str
		ENDELSE
	END
    ELSE:
    ENDCASE

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcemrirs_Redraw_im1, Pst

WIDGET_CONTROL, (*Pst).wd.wDraw1,  GET_VALUE = window_draw
window_draw -> Draw, (*Pst).vi1.oView

RETURN, 1

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwdcmri_GetResultBaseName, Pst, PARAMETER_IMAGE=parameter_image

	model     = (*(*Pst).modelinfo).model
	parameter = (*(*Pst).modelinfo).ids[(*Pst).opt_current_param]

	IF ((*Pst).flag_stddev EQ 1) THEN str_stddev = '.stddev' ELSE str_stddev = ''

	str_pointini = fwdcemri_suffixresultname(ST_ROI=(*Pst).roi, SIZEIMAGE=(*Pst).inf.size_im, '[Im]')

	IF KEYWORD_SET(parameter_image) THEN $
		str_file_ini = 'Param_' + model + '.' + parameter + str_stddev + '.' + str_pointini

	RETURN, str_file_ini

END

;**************************************************************************************************
;**************************************************************************************************

PRO fwdcemrirs_ViewParametricMap, Pst

	n_points_lr = N_ELEMENTS(*(*Pst).roi.idx_points_lr)
   	data_Signal_lr = (*(*Pst).results.dataptr[(*Pst).opt_current_param,(*Pst).flag_stddev])[0:n_points_lr-1]
	range  = REFORM((*Pst).results.scales[*,(*Pst).opt_current_param, (*Pst).flag_stddev])
	;-----------------------------------------------------------------------------------------
	str_colorbar_title = (*(*Pst).modelinfo).idsinf[2,(*Pst).opt_current_param]
	IF (*Pst).flag_stddev EQ 1 THEN str_colorbar_title+=' (SD) '
    ;-----------------------------------------------------------------------------------------

	IF N_ELEMENTS(data_Signal_lr) GE 1 THEN BEGIN

		(*Pst).vi1.opaletteROI->getproperty, RED_VALUES=r, GREEN_VALUES=g, BLUE_VALUES=b

		minval_lr = range[0]
		maxval_lr = range[1]
		IF maxval_lr LE minval_lr THEN BEGIN
			maxval_lr+=ABS(maxval_lr)/100.0 ; maxval !=minval
			minval_lr-=ABS(minval_lr)/100.0 ; maxval !=minval
		ENDIF
		IF maxval_lr EQ minval_lr THEN BEGIN ; si los dos eran cero
			maxval_lr+=1
			minval_lr-=1
		ENDIF
		WIDGET_CONTROL, (*Pst).wd.wtext_palette_minval, SET_VALUE=STRTRIM(STRING(minval_lr, FORMAT='(8f0)'),2)
		WIDGET_CONTROL, (*Pst).wd.wtext_palette_maxval, SET_VALUE=STRTRIM(STRING(maxval_lr, FORMAT='(8f0)'),2)
		;---------------------------------------------------------
		data_lr_scl = BYTSCL(data_Signal_lr[0:N_ELEMENTS(*(*Pst).roi.idx_points_lr)-1], MIN=minval_lr, MAX=maxval_lr)
		;---------------------------------------------------------

		CASE (*Pst).opt_drawROIoutliers OF
		1 : BEGIN ; transparent
			pos_nodraw = WHERE((data_Signal_lr GT maxval_lr) OR  (data_Signal_lr LT minval_lr),ct_nodraw)
			IF (ct_nodraw EQ 0) THEN undefine, pos_nodraw
		END
		2 : BEGIN ; black
			pos_drawblack = WHERE((data_Signal_lr GT maxval_lr) OR  (data_Signal_lr LT minval_lr),ct_drawblack)
			IF (ct_drawblack EQ 0) THEN undefine, pos_drawblack
		END
		ELSE : ; 0
		ENDCASE

		imag_lr = imagFromArray(data_lr_scl, DIMENSIONS=(*Pst).roi.roi_sz/(*Pst).roi.np_lr,$
			 IDX_POINTS = *(*Pst).roi.idx_points_lr, PINI=(*Pst).roi.pr_ini, MASK=mask, $
			 IDX_NOVALID= pos_nodraw)
		img_color = imag2threeChannels(imag_lr, RED_VALUES=r, GREEN_VALUES=g, BLUE_VALUES=b)

		IF (N_ELEMENTS(pos_drawblack) GT 0) THEN BEGIN
			pos2_drawblack = (*(*Pst).roi.idx_points_lr)[pos_drawblack]
			size_xy = (SIZE(img_color, /DIMENSIONS))[0]*((SIZE(img_color, /DIMENSIONS))[1])
			img_color[pos2_drawblack]=0
			img_color[pos2_drawblack+size_xy]=0
			img_color[pos2_drawblack+size_xy*2]=0
		ENDIF
		;---------------------------------------------------------
		img_alphachannel = [[[img_color]],[[mask*255B]]]
		dimensions = SIZE(mask, /DIMENSIONS)
		(*Pst).vi1.oimageROI -> setProperty, DATA = img_alphachannel, DIMENSIONS=(*Pst).roi.roi_sz, $
			INTERLEAVE=2, GREYSCALE=0, LOCATION=(*Pst).roi.pr_ini, HIDE=0

		(*Pst).vi1.ohcolorbar->setproperty, RANGE=[minval_lr,maxval_lr], PALETTE=(*Pst).vi1.opaletteROI,$
    		TITLE=str_colorbar_title

		(*Pst).vi1.opolyline_ROI->setproperty, HIDE = (*Pst).flag_roilimits EQ 0
		WIDGET_CONTROL, (*Pst).wd.wDraw1, GET_VALUE = window_draw
    	window_draw -> Draw, (*Pst).vi1.oView

	ENDIF ELSE BEGIN
		(*Pst).vi1.ohcolorbar-> setproperty, PALETTE=(*Pst).vi1.opalette, RANGE=[0,255], TITLE=''
		(*Pst).vi1.oimageROI -> setProperty, HIDE=1
		WIDGET_CONTROL, (*Pst).wd.wDraw1, GET_VALUE = window_draw
		(*Pst).vi1.opolyline_ROI->setproperty, HIDE = 0
    	window_draw -> Draw, (*Pst).vi1.oView
	ENDELSE

END


;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwBaseStruct_Transparency, PARENT_BASE=parent_base, MAIN_BASE=main_base, FRAME=frame


main_base = WIDGET_BASE(parent_base, /COLUMN, /ALIGN_CENTER, /BASE_ALIGN_CENTER, FRAME=frame)
	base_1a = WIDGET_BASE(main_base, /ROW, /ALIGN_LEFT, /BASE_ALIGN_CENTER)
		base_1a1 = WIDGET_BASE(base_1a, /ROW, /ALIGN_CENTER, /BASE_ALIGN_RIGHT)
		base_1a2 = WIDGET_BASE(base_1a, /ROW, /ALIGN_CENTER, /BASE_ALIGN_RIGHT)
		base_1a3 = WIDGET_BASE(base_1a, /ROW, /ALIGN_CENTER, /BASE_ALIGN_RIGHT, /NONEXCLUSIVE)
	base_1b = WIDGET_BASE(main_base, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
		base_1b1 = WIDGET_BASE(base_1b , /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_1b2 = WIDGET_BASE(base_1b , /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_1b3 = WIDGET_BASE(base_1b , /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)


strarr_drawroioutliers = ['Palette', 'Transp.', 'Black']
dropl_roidrawmode = WIDGET_DROPLIST(base_1a1, TITLE='Outliers', UNAME='DROPLIST_DRAWROIOUTLIERS',$
	UVALUE='DROPLIST_DRAWROIOUTLIERS', VALUE=strarr_drawroioutliers, XSIZE=110)

slide_transparency = WIDGET_SLIDER(base_1a2, MAXIMUM=100, MINIMUM=0, SCR_XSIZE=size_slider_3, $
    VALUE=0, UVALUE='TRANSPARENCY_CHANGE', UNAME = 'SLICE_TRANSPARENCY',  $
    TITLE='Transparency(%)', /SUPPRESS_VALUE)
label_transparency  = WIDGET_LABEL(base_1a2, VALUE='0%', SCR_XSIZE=25, SCR_YSIZE=20)


label_palette_minval = WIDGET_LABEL(base_1b1, VALUE='Scale: Min')
text_palette_minval  = WIDGET_TEXT(base_1b1,  VALUE= '',  UVALUE='ADJUST_PALETTE', $
	/EDITABLE, XSIZE=8, UNAME= 'PALETTE_MINVAL')
label_palette_maxval = WIDGET_LABEL(base_1b2, VALUE='Max')
text_palette_maxval  = WIDGET_TEXT(base_1b2,  VALUE= '',  UVALUE='ADJUST_PALETTE', $
	/EDITABLE, XSIZE=8, UNAME= 'PALETTE_MAXVAL')
bttn_adjust = WIDGET_BUTTON(base_1b3, VALUE='Reset', UVALUE='ADJUST_PALETTE', UNAME='RESET')

st_widget_transparency = { $
		dropl_roidrawmode : dropl_roidrawmode,  $
		slide_transparency:	slide_transparency, $
		label_transparency: label_transparency, $
		text_palette_minval:text_palette_minval,$
		text_palette_maxval:text_palette_maxval,$
		bttn_adjust : bttn_adjust }

RETURN, st_widget_transparency

END


;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwBaseStruct_KineticParameter, PARENT_BASE=parent_base, MAIN_BASE=main_base, FRAME=frame, TITLE=title, VALUES=values

main_base = WIDGET_BASE(parent_base, /ROW,/ALIGN_CENTER, /BASE_ALIGN_CENTER, FRAME=frame)
	base_1 = WIDGET_BASE(main_base, /ROW, /ALIGN_LEFT,  /BASE_ALIGN_LEFT)
	base_2 = WIDGET_BASE(main_base, /ROW, /ALIGN_RIGHT, /BASE_ALIGN_RIGHT, /EXCLUSIVE)

dropls_params = WIDGET_DROPLIST(base_1, TITLE=title,   UNAME= 'DROPLIST_MODEL',UVALUE='VIEW_PARAM', VALUE=values)
bttn_value    = WIDGET_BUTTON(base_2, VALUE= 'Value',  UNAME='BUTTON_ACTIVATE_VALUE', UVALUE='VIEW_PARAM')
bttn_stddev   = WIDGET_BUTTON(base_2, VALUE= 'Std.dev',UNAME='BUTTON_ACTIVATE_STDDEV',UVALUE='VIEW_PARAM')

st_widgets_parameter = {$
					dropls_params:dropls_params,$
					bttn_value  : bttn_value, $
					bttn_stddev : bttn_stddev}

RETURN, st_widgets_parameter

END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_DCEMRIresult_destroy, ev


	WIDGET_CONTROL, ev, GET_UVALUE=pst

	PTR_FREE, (*Pst).data
	PTR_FREE, (*Pst).mask
	PTR_FREE, (*Pst).parametric_AIF
	PTR_FREE, (*Pst).results.arr_Cp
	PTR_FREE, (*Pst).results.arr_Signal_AIF
	PTR_FREE, (*Pst).results.arr_signal_RR
	PTR_FREE, (*Pst).results.arr_time
	PTR_FREE, (*Pst).results.data_CT
	PTR_FREE, (*Pst).results.arr_Ct_RR
	PTR_FREE, (*Pst).results.data_Signal_lr
	PTR_FREE, (*Pst).results.data_Signal_ROI
	PTR_FREE, (*Pst).roi.mask_lr
	PTR_FREE, (*Pst).roi.idx_points
	PTR_FREE, (*Pst).roi.idx_points_lr
	PTR_FREE, (*Pst).roi.im_idxs

	IF OBJ_VALID((*Pst).vi1.ocontainer)  THEN BEGIN
	    OBJ_DESTROY, (*Pst).vi1.ocontainer
	ENDIF

	HEAP_GC
	WIDGET_CONTROL,  (*Pst).wd.wbase_main,  /DESTROY
	HEAP_GC

END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_DCEMRIresult_event, ev

opt_catcherror = 0
;*****************************************
IF LMGR(/RUNTIME) OR LMGR(/VM) OR (opt_catcherror EQ 1) THEN BEGIN
	CATCH, theError
	IF theError NE 0 THEN BEGIN
	  CATCH, /Cancel
	  void = DIALOG_MESSAGE(!ERROR_STATE.MSG, /ERROR)
	  RETURN
	ENDIF
ENDIF
;*****************************************

WIDGET_CONTROL, ev.top, GET_UVALUE=Pst
WIDGET_CONTROL, ev.id,  GET_UVALUE=uval
HEAP_GC

IF N_ELEMENTS(uval) EQ 0 THEN RETURN

CASE STRUPCASE(uval[0]) OF
;-------------------------------------------------------------------------------------------------------
'BABOUT' : BEGIN

	ok = Interface_DCEMRI_About(GROUP_LEADER=(*Pst).wd.wbase_main)
	RETURN
END
;-------------------------------------------------------------------------------------------------------
'BQUIT' : BEGIN

    Interface_DCEMRIresult_destroy, ev.top
    HEAP_GC &  RETURN
END ;BQUIT
;-------------------------------------------------------------------------------------------------------
'BUTTON_ACTIVATE_STDDEV': BEGIN

	(*Pst).flag_stddev = WIDGET_INFO((*Pst).wd.wbttn_stddev, /BUTTON_SET)

	number = WIDGET_INFO((*Pst).wd.dropls_params, /DROPLIST_SELECT)
	WIDGET_CONTROL, (*Pst).wd.dropls_params, SEND_EVENT={ID:0L, TOP:0L, HANDLER:0L}
END
'CLOSE_WINDOWS' : BEGIN
	WDEL_ALL
	RETURN
END
;-------------------------------------------------------------------------------------------------------
'DROPLIST_DRAWROIOUTLIERS' : BEGIN

	(*Pst).opt_drawROIoutliers = WIDGET_INFO((*Pst).wd.wdropl_roidrawmode, /DROPLIST_SELECT)
	WIDGET_CONTROL, (*Pst).wd.dropls_params, SEND_EVENT={ID:0L, TOP:0L, HANDLER:0L}
END
;-------------------------------------------------------------------------------------------------------
'INTERPOLATION': BEGIN

     (*Pst).flag_interp = WIDGET_INFO((*Pst).wd.wbttn_interp, /BUTTON_SET)
     ok=fwdcemri_Draw_Image(Pst, NODRAW=0)
END
;-------------------------------------------------------------------------------------------------------
'ROI_LIMITS' : BEGIN

	(*Pst).flag_roilimits = WIDGET_INFO((*Pst).wd.bttn_roilimits, /BUTTON_SET)
	(*Pst).vi1.opolyline_ROI->setproperty, HIDE = (*Pst).flag_roilimits EQ 0
    ok = fwdcemrirs_Redraw_im1(Pst)
END
;-------------------------------------------------------------------------------------------------
'S_CHANGE' : BEGIN   ; Elige la imagen 2D a mostrar (cambiando slice o frame)

    uname= WIDGET_INFO(ev.id, /UNAME)
    ok = fwdcemri_Draw_Image(Pst)
END
;-------------------------------------------------------------------------------------------------------
'PLOTSOPENED_TYPE' : BEGIN

	(*Pst).opt_plotsOpened = fwdcemri_PlotsOpenedType(ev.id, (*Pst).wd.wbttn_plotsOpened, (*Pst).opt_plotsOpened)

END
;--------------------------------------------------------------------------------------------------
'VIEW_PARAM' : BEGIN

	IF (ev.id EQ (*Pst).wd.dropls_params) THEN BEGIN ; identificador de droplist activado
		number_param = WIDGET_INFO(ev.id, /DROPLIST_SELECT)
		WIDGET_CONTROL, ev.id, SET_DROPLIST_SELECT=number_param, /UPDATE
		str_paramidentifier = (*(*Pst).modelinfo).ids[number_param]
		(*Pst).opt_current_param = number_param

		opt_stddev     = WIDGET_INFO((*Pst).wd.wbttn_stddev, /BUTTON_SET)

		(*Pst).flag_stddev = (*(*Pst).modelinfo).sd_valid[(*Pst).opt_current_param] AND opt_stddev
		WIDGET_CONTROL, (*Pst).wd.wbttn_stddev, SET_BUTTON=(*Pst).flag_stddev
		WIDGET_CONTROL, (*Pst).wd.wbttn_stddev, SENSITIVE= (*(*Pst).modelinfo).sd_valid[(*Pst).opt_current_param]

	ENDIF

	IF (ev.id EQ (*Pst).wd.wbttn_stddev) THEN BEGIN ; identificador de droplist activado
		opt_stddev     = WIDGET_INFO((*Pst).wd.wbttn_stddev, /BUTTON_SET)
		(*Pst).flag_stddev = (*(*Pst).modelinfo).sd_valid[(*Pst).opt_current_param] AND opt_stddev
		WIDGET_CONTROL, (*Pst).wd.wbttn_value, SET_BUTTON= (*Pst).flag_stddev EQ 0
	ENDIF

	IF N_ELEMENTS(*(*Pst).results.dataptr[(*Pst).opt_current_param,(*Pst).flag_stddev]) LE 0 THEN BEGIN
		ok = DIALOG_MESSAGE('Value not calculated', /INFO)
		RETURN
	ENDIF

	fwdcemrirs_ViewParametricMap, Pst

END
;-------------------------------------------------------------------------------------------------------
'TRANSPARENCY_CHANGE' : BEGIN
	WIDGET_CONTROL, (*Pst).wd.wslide_transparency, GET_VALUE=transparency_pc
	(*Pst).vi1.oimageROI->setproperty, ALPHA_CHANNEL=1-transparency_pc/100.0;, BLEND_FUNCTION=[1,4]
	WIDGET_CONTROL,(*Pst).wd.wlabel_transparency, SET_VALUE=STRTRIM(STRING(transparency_pc),2) + '%'
	ok=fwdcemri_Draw_Image(Pst)
END
;-------------------------------------------------------------------------------------------------------
'PALETTE1' : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)
	(*Pst).vi1.opalette   -> Getproperty, RED=r,GREEN=G, BLUE=b
	TVLCT, r,g,b
 	XLOADCT, GROUP=(*Pst).wd.wbase_main, /MODAL, UPDATECALLBACK='fwdcemrirs_Change_palette1b', UPDATECBDATA=*Pst, /USE_CURRENT
 END
;-------------------------------------------------------------------------------------------------------
'PALETTE3' : BEGIN

	(*Pst).vi1.opaletteROI   -> Getproperty, RED=r,GREEN=G, BLUE=b
	TVLCT, r,g,b
	XLOADCT, GROUP=(*Pst).wd.wbase_main, /MODAL, UPDATECALLBACK='fwdcemrirs_Change_palette2b', UPDATECBDATA=*Pst, /USE_CURRENT
	fwdcemrirs_ViewParametricMap, Pst
END
;-------------------------------------------------------------------------------------------------------
'ADJUST_PALETTE' : BEGIN

	identifier = (*(*Pst).modelinfo).ids[(*Pst).opt_current_param]
	iddata = (*Pst).opt_current_param

	uname= WIDGET_INFO(ev.id, /UNAME)
    CASE uname OF
    'RESET' : BEGIN
    	(*Pst).results.scales[*,iddata,(*Pst).flag_stddev] = [MIN(*(*Pst).results.dataptr[iddata,(*Pst).flag_stddev]), $
    		MAX(*(*Pst).results.dataptr[iddata, (*Pst).flag_stddev])]
	END
	ELSE: BEGIN

		WIDGET_CONTROL, (*Pst).wd.wtext_palette_minval, GET_VALUE=min_value
		WIDGET_CONTROL, (*Pst).wd.wtext_palette_maxval, GET_VALUE=max_value
		fmin_value = (FLOAT(min_value))[0]
		fmax_value = (FLOAT(max_value))[0]
		IF fmin_value GE fmax_value THEN BEGIN
			ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main,'Error, min value must be less then max value, ', /INFO)
			RETURN
		ENDIF
		WIDGET_CONTROL, (*Pst).wd.wtext_palette_minval, SET_VALUE=STRTRIM(STRING(fmin_value, FORMAT='(8f0)'),2)
		WIDGET_CONTROL, (*Pst).wd.wtext_palette_maxval, SET_VALUE=STRTRIM(STRING(fmax_value, FORMAT='(8f0)'),2)

		(*Pst).results.scales[*,iddata, (*Pst).flag_stddev] = [fmin_value,fmax_value]

	END
	ENDCASE

	WIDGET_CONTROL, (*Pst).wd.dropls_params, SEND_EVENT={ID:0L, TOP:0L, HANDLER:0L}
END
;-------------------------------------------------------------------------------------------------------
'BACKGROUND_IMAGE_SIGNAL' : BEGIN

	WIDGET_CONTROL, (*Pst).wd.wDraw1, GET_VALUE = window_draw
	CASE (*Pst).opt.background_MRI OF
	1 : BEGIN ; first color of palette
		(*Pst).vi1.opalette->getproperty, RED_VALUES=r, GREEN_VALUES=g, BLUE_VALUES=b
		(*Pst).vi1.oview->setproperty, COLOR=[r[0],g[0],b[0]]
		(*Pst).vi1.ohcolorbar->setproperty, COLOR= (*Pst).opt.colorbar_labelcolor
		(*Pst).opt.background_MRI = 2
	END
	2 : BEGIN ; default color
		(*Pst).vi1.oview->setproperty, COLOR=[255,255,255L]
		(*Pst).vi1.ohcolorbar->setproperty, COLOR= [0,0,0d]
		(*Pst).opt.background_MRI = 3
	END
	3 : BEGIN ; white color (and numbers of palette in black)
		(*Pst).vi1.oview->setproperty, COLOR=[240,240,240L]
		(*Pst).vi1.ohcolorbar->setproperty, COLOR= [0,0,0d]
		(*Pst).opt.background_MRI = 4
	END
	4 : BEGIN ; default color
		(*Pst).vi1.oview->setproperty, COLOR=(*Pst).opt.background_color_mri
		(*Pst).vi1.ohcolorbar->setproperty, COLOR= (*Pst).opt.colorbar_labelcolor
		(*Pst).opt.background_MRI = 1
	END
	ENDCASE

    window_draw -> Draw, (*Pst).vi1.oView
END
;-------------------------------------------------------------------------------------------------------
'COLORBAR_FONT' : BEGIN

	uname= WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'INCREASE': val =  1
	'DECREASE': val = -1
	ENDCASE

	(*Pst).vi1.ohColorbar->Getproperty, FONTSIZE= fontsize
	(*Pst).vi1.ohColorbar->Setproperty, FONTSIZE= 1 > (fontsize+val)

	WIDGET_CONTROL, (*Pst).wd.wDraw1, GET_VALUE = window_draw
	window_draw -> Draw, (*Pst).vi1.oView
END
;-------------------------------------------------------------------------------------------------------

'ROI_MODEL' : BEGIN

	n_hd = 1l
	;model_id = (*(*Pst).modelinfo).idsinf[4,(*Pst).opt_current_param]

	IF (*(*Pst).modelinfo).model EQ 'CURVE' THEN BEGIN
		ok = DIALOG_MESSAGE('Semi-parametric parameters: No kinetic model', /INFO, DIALOG_PARENT=ev.top) & RETURN
	ENDIF
	p_idx  = N_ELEMENTS(*(*Pst).roi.idx_points_lr); last element is the ROI
	IF p_idx LE 0 THEN BEGIN
		ok = DIALOG_MESSAGE('ROI analysis must be done before', /INFO, DIALOG_PARENT=ev.top) & RETURN
	ENDIF

	nparam       = (*(*Pst).modelinfo).nparam
	parameters   = DBLARR(nparam)
	standard_dev = DBLARR(nparam)

	uname= WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'VIEW_ROI_MODEL' :BEGIN
		nparam = (*(*Pst).modelinfo).nparam
		parameters   = DBLARR(nparam)
		standard_dev = DBLARR(nparam)
		FOR i=0l, nparam-1 DO BEGIN
			parameters[i]     =  (*(*Pst).results.dataptr[i,0])[p_idx]
			IF (*(*Pst).modelinfo).sd_valid[i] THEN $
				standard_dev[i]   =  (*(*Pst).results.dataptr[i,1])[p_idx]
		ENDFOR
	END
	;'VIEW_ROI_MODEL_NOOUTLIERS' : BEGIN
	;	idx_points = fwdcemri_RemakeIndexWithoutOutliers(Pst, MASK_LR=mask_lr)
	;	arr_data   = fwdcemri_analyze_ROIwithOutliers(Pst, MODEL_ID=model_st, IDX_POINTS=idx_points,  PARAMS=params, STANDARD_DEV=standard_dev, ARR_CT=arr_ct)
	;	*(*Pst).results.data_ROI_outliers = arr_data
	;	parameters = params
	;	undefine, p_idx
	;END
	ENDCASE

	win_number  = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)

	ok = fwdcemrirs_Plot_Model(Pst,  WIN_NUMBER=win_number, POINT=point, $
		P_IDX=p_idx, XPOS=xpos, YPOS= ypos, N_HD=n_hd, $
		ARR_DATA=arr_data, PARAMETERS = parameters, STANDARD_DEV=standard_dev)

END
;--------------------------------------------------------------------------------------------------
'EXPORT_IMAGE' : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)
	CASE uname OF
	'EXPORT_WINDOW' : BEGIN
		WIDGET_CONTROL, (*Pst).wd.wDraw1,  GET_VALUE = window_1
		window_1->getproperty, IMAGE_DATA=image_view
	END
	'EXPORT_IMAGE'    : BEGIN
		(*Pst).vi1.oimage -> GetProperty, DATA = image, DIMENSIONS=dimensions, LOCATION=location
		;(*Pst).vi1.hcolorbar -> GetProperty, DATA = image, DIMENSIONS=dimensions, LOCATION=location
	END
	'EXPORT_IMAGEROI' : BEGIN
		scl=6

		(*Pst).vi1.oimageROI -> GetProperty, DATA = img_alphachannel, DIMENSIONS=dimensions, LOCATION=location
		roi_color  = TEMPORARY(img_alphachannel[*,*,0:2])
		image_view = REBIN(roi_color, dimensions[0]*scl,dimensions[1]*scl,3, /SAMPLE)
		image_view = change_true(image_view)
		opt_no_view = 1
	END
	ENDCASE

	str_file_ini = fwdcmri_GetResultBaseName(Pst, /PARAMETER_IMAGE) + '.ROI'

	path_result = (*Pst).path_result
	im = fwdcemri_exportImage(image_view, NAME=str_file_ini, OPTION_NAME=0, INFO=*(*Pst).inf.info, PATH_RESULT=path_result, NO_VIEW=opt_no_view)
	(*Pst).path_result = path_result

END
;-------------------------------------------------------------------------------------------------------
'ROI_KINETICS' : BEGIN

	uname= WIDGET_INFO(ev.id, /UNAME)
	flag_errorbars = 1
	IF flag_errorbars THEN type_plot = 2 ELSE type_plot = 1
	CASE uname OF
	'VIEW_ROI_KINETICS': BEGIN
		win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
		data_kinetics = fwdcemri_ROIkinetics(Pst, PLOT=1, WIN_NUMBER=win_number, XPOS=xpos, YPOS=ypos, TYPE_PLOT=type_plot)
	END
	'VIEW_ROI_KINETICS_NOOUTLIERS' : BEGIN
		;IF (*Pst).flag_analysis_done EQ 0 THEN BEGIN
		;	ok = DIALOG_MESSAGE('ROI analysis must be done', /INFO) & RETURN
		;ENDIF
		idx_points    = fwdcemri_RemakeIndexWithoutOutliers(Pst, MASK_LR=mask_lr)
		win_number    = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
		data_kinetics = fwdcemri_ROIkinetics(Pst, IDX_POINTS=idx_points, PLOT=1, WIN_NUMBER=win_number, XPOS=xpos, YPOS=ypos, TYPE_PLOT=type_plot)
	END
	ELSE:
	ENDCASE
END
;-------------------------------------------------------------------------------------------------------
'EXPORT_ASCII_ROI_MULTICOLUMN' : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)

	reverse_lines = 1 ; por defecto
	str_info = ''

	n_points     = N_ELEMENTS(*(*Pst).roi.idx_points_lr)
	n_parameters =  (*(*Pst).modelinfo).nparam

	str_info  = STRARR(60)
	str_info[0] = '## Model: ' + (*(*Pst).modelinfo).str_model

	n=2l
	FOR i=0l, n_parameters-1 DO BEGIN
		im_valtemp          = (*(*Pst).results.dataptr[i,0])[0:n_points-1]
		min_val_allowed     = (*Pst).results.scales[0,i,0]
		max_val_allowed     = (*Pst).results.scales[1,i,0]
		str_info[n++]       = '## Column ' + STRTRIM(n-2,2)  + ': ' + (*(*Pst).modelinfo).idsinf[2,i] + ' - ' + (*(*Pst).modelinfo).idsinf[3,i]

		IF (*(*Pst).modelinfo).sd_valid[i] THEN BEGIN
			im_stddev           = (*(*Pst).results.dataptr[i,1])[0:n_points-1]
			max_stddev_allowed  = (*Pst).results.scales[1,i,1]

			str_info[n++]   = '## Column ' + STRTRIM(n-2,2) + ': '+  (*(*Pst).modelinfo).idsinf[2,i] + ' (Estimated SD) '

			IF i EQ 0 THEN BEGIN
				im_valtotal=[[im_valtemp],[im_stddev]]
			ENDIF ELSE BEGIN
				im_valtotal=[[im_valtotal],[im_valtemp],[im_stddev]]
			ENDELSE
		ENDIF ELSE BEGIN
			IF i EQ 0 THEN BEGIN
				im_valtotal=im_valtemp
			ENDIF ELSE BEGIN
				im_valtotal=[[im_valtotal],[im_valtemp]]
			ENDELSE
		ENDELSE
	ENDFOR

	str_info = str_info[0:n+6]
	str_info[n++] = '## Column ' + STRTRIM(n-2,2) + ': ' +  'Pixel index (x-y)'
	n+=2
	str_info[n++] = '## Slice size (x,y): ' + STRTRIM((*Pst).inf.size_im[0],2) +  ' ' + STRTRIM((*Pst).inf.size_im[1],2)
	str_info[n++] = '## ROI resolution  : ' + STRTRIM((*Pst).roi.np_lr[0],2)

	array_multcolumn = ROTATE(im_valtotal,4)

	IF n_points EQ 1 THEN array_multcolumn = REFORM(array_multcolumn, N_ELEMENTS(array_multcolumn),1)

	str_array = Translate_data2strings(array_multcolumn, N_DECIMALS=5, REVERSE_LINES=reverse_lines, $
		MIN_COL_LENGTH=15, INTEGER=opt_integer, MAX_VALUE=arr_maxvalues, STR_NOVALID='N.D',	COL_TITLES=col_titles)

	IF (SIZE(str_array,/TYPE) NE 7) THEN BEGIN
		IF str_array[0] EQ -1 THEN BEGIN
			ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main, 'Nothing to export', /INFO) &	RETURN
		ENDIF
	ENDIF

	idx_points_corners = IdxPoints_FirstFromlowRes(MASK_IDX=*(*Pst).roi.im_idxs, DIMENSIONS=(*Pst).inf.size_im, $
			ROI_RESOLUTION=(*Pst).roi.np_lr, MASK_OUTPUT=mask_output)
	idx_points_Yup = IdxPoints_ChangeYconvention(idx_points_corners, DIMENSIONS=(*Pst).inf.size_im, ARRAY_XY=array_xy)

	IF n_points EQ 1 THEN array_xy = REFORM(array_xy, N_ELEMENTS(array_xy),1)
	str_idx_points    = Translate_data2strings(array_xy, N_DECIMALS=5, REVERSE_LINES=reverse_lines, INTEGER=1)
	;str_idx_points    = ['Index (x-y)',REFORM(str_idx_points)]
	FOR i=0l, N_ELEMENTS(str_array)-1 DO BEGIN
		str_array[i] = str_array[i] + '    ' + str_idx_points[i]
	ENDFOR

	str_array = [str_info, str_array]

	IF N_ELEMENTS(str_label) NE 0 THEN str_array = [str_label, str_array]

	str_pointini = fwdcemri_suffixresultname(ST_ROI=(*Pst).roi, SIZEIMAGE=(*Pst).inf.size_im, '[ROImc]')

	str_file_ini = 'Model_' + STRLOWCASE((*(*Pst).modelinfo).model) + '_MultiColumn_'
	str_file_ini+=str_pointini +'.txt

	str_file = DIALOG_PICKFILE(PATH=(*Pst).path_result, TITLE='Export model params (ASCII file)', $
		DEFAULT_EXTENSION = 'txt', WRITE=1, FILE=str_file_ini, DIALOG_PARENT=ev.top)
	IF str_file EQ '' THEN RETURN

	(*Pst).path_result = get_path(str_file)
	file_name  = get_path(str_file) + get_name_field(str_file) + '.txt'
	ok = WriteFile_ASCII(file_name, ARRAY=str_array, ERROR_ID=error_id)

	IF ok EQ 1 THEN BEGIN
		ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main,['Result saved in file ', get_name(file_name)], /INFO)
	ENDIF
END
;-------------------------------------------------------------------------------------------------------
'EXPORT_ASCII' : BEGIN

	uname = WIDGET_INFO(ev.id, /UNAME)

	reverse_lines = 1     ; por defecto
	opt_index_points = 1L ; opción de poner los indices de los puntos
	opt_name         = 0
	str_info = ''

	CASE uname OF
	'EXPORT_ROI_RESULTS' : BEGIN
		identifier = (*(*Pst).modelinfo).ids[(*Pst).opt_current_param]
		iddata     = (*Pst).opt_current_param

		n_points = N_ELEMENTS(*(*Pst).roi.idx_points_lr)
		i_stddev = (*Pst).flag_stddev*(*(*Pst).modelinfo).sd_valid[iddata] ;se escribe o no la desviación estándar
		;------------------------------------
		im_temp  = (*(*Pst).results.dataptr[iddata,i_stddev])[0:n_points-1]
		str_info = (*(*Pst).modelinfo).idsinf[2,iddata]
		;------------------------------------
		str_file_ini = 'Param_' + identifier + '_'
		IF i_stddev EQ 1 THEN str_file_ini+='SD_'
		;------------------------------------
		IF uval[1] EQ 'MATRIX_FORMAT' THEN BEGIN
    		IF ((*Pst).roi.opt_typeROI NE 0) THEN BEGIN
    			ok = DIALOG_MESSAGE('Matrix format output not allowed with this ROI type', /INFO)
    			RETURN
    		ENDIF
    		im_matrix = FLTARR((*Pst).inf.size_roi/(*Pst).roi.np_lr)
    		im_matrix[*(*Pst).roi.idx_points_lr] = im_temp
    		im_temp = TEMPORARY(im_matrix)
    		opt_name      = 1
    		str_file_ini+='MatrixF_'
    	ENDIF
    END
    'EXPORT_KINETIC_CURVE' : BEGIN  ; cinética temporal
    	win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
		im_temp = fwdcemri_ROIkinetics(Pst, WIN_NUMBER=win_number,XPOS=xpos, YPOS=ypos);
		str_file_ini = 'ROI_kinetics'
    	reverse_lines = 0 ;
    	opt_index_points = 0l
		opt_name      = 3
		str_label = '  Time (s) ' + STRING(9B) + 'Average (AU) ' + STRING(9B) + 'SD'
		undefine, identifier
	END
	'EXPORT_KINETIC_CURVE_OUTLIERS' : BEGIN
		idx_points = fwdcemri_RemakeIndexWithoutOutliers(Pst, MASK_LR=mask_lr)
		win_number = get_windowNumber(OPTNEXT=(*Pst).opt_plotsOpened EQ 1, XPOS=xpos, YPOS=ypos)
		im_temp = fwdcemri_ROIkinetics(Pst, IDX_POINTS=idx_points, WIN_NUMBER=win_number,XPOS=xpos, YPOS=ypos)
		str_file_ini = 'ROI_kinetics_NoOutliers'
    	reverse_lines = 0 ;
    	opt_index_points = 0l
		opt_name      = 3
		str_label = '  Time (s) ' + STRING(9B) + 'Average (AU) ' + STRING(9B) + 'SD'
		undefine, identifier
	END
	ELSE:
    ENDCASE
    ;----------------------------------------------------------
	str_array = Translate_data2strings(im_temp, N_DECIMALS=5, REVERSE_LINES=reverse_lines, INTEGER=opt_integer)

	IF (SIZE(str_array,/TYPE) NE 7) THEN BEGIN
		IF str_array[0] EQ -1 THEN BEGIN
			ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main, 'Nothing to export', /INFO)
			RETURN
		ENDIF
	ENDIF
	IF N_ELEMENTS(str_label) NE 0 THEN str_array = [str_label, str_array]

	str_pointini = fwdcemri_suffixresultname(ST_ROI=(*Pst).roi, SIZEIMAGE=(*Pst).inf.size_im, '[ROIf]')

	str_file_ini+=str_pointini +'.txt

	str_file = DIALOG_PICKFILE(PATH=(*Pst).path_result, TITLE='Export ' + str_info  + ' (ASCII file)', $
		DEFAULT_EXTENSION = 'txt', WRITE=1, FILE=str_file_ini, DIALOG_PARENT=ev.top)
	IF str_file EQ '' THEN RETURN

	(*Pst).path_result = get_path(str_file)
	file_name  = get_path(str_file) + get_name_field(str_file) + '.txt'
	ok = WriteFile_ASCII(file_name, ARRAY=str_array, ERROR_ID=error_id)

	IF ok EQ 1 THEN BEGIN
		ok = DIALOG_MESSAGE(DIALOG_PARENT=(*Pst).wd.wbase_main,['Result saved in file ', get_name(file_name)], /INFO)
	ENDIF

END
;-------------------------------------------------------------------------------------------------------
'EXPORT_ROI' : BEGIN

	path_results = (*Pst).path_result
	ok = fwdcemri_ExportROI((*Pst).roi, IMAGE_SIZE = (*Pst).inf.size_im, PATH=path_results , OPOLYLINE= (*Pst).vi1.opolyline_ROI)
	(*Pst).path_result =path_results
END
;-------------------------------------------------------------------------------------------------------
'DRAW' : BEGIN     ; Drawing events

	ok = fwdcemrirs_DrawingEvents_Result(ev, Pst)
END
;-------------------------------------------------------------------------------------------------------
'ZOOM_CHANGE' : BEGIN     ; Cambia el zoom

    uname= WIDGET_INFO(ev.id, /UNAME)
    CASE uname OF
    'ZOOM_MENOS' : BEGIN
    	ok = fwdcemri_Draw_Image(Pst, CONTROL_ZOOM=0.8d,  NODRAW=1)
    END
    'ZOOM_MAS'   : BEGIN
    	ok = fwdcemri_Draw_Image(Pst, CONTROL_ZOOM=1.25d, NODRAW=1)
    END
    ENDCASE
    ok = fwdcemrirs_Redraw_im1(Pst)
END
;-------------------------------------------------------------------------------------------------------
ELSE : ;PRINT, '-'
ENDCASE

;HELP, /MEMORY

END

;**************************************************************************************************
;**************************************************************************************************

PRO  Interface_DCEMRIresult, Pst, TITLE=title, WIDGET_ID=widget_id

COMPILE_OPT DEFINT32, STRICTARR, HIDDEN

opt_catcherror = 0
;*****************************************
IF LMGR(/RUNTIME) OR LMGR(/VM) OR (opt_catcherror EQ 1) THEN BEGIN
	Catch, theError
	IF theError NE 0 THEN BEGIN
		Catch, /Cancel
		void = DIALOG_MESSAGE(!ERROR_STATE.MSG, /ERROR)
		;void = ERROR_MESSAGE()
		RETURN
	ENDIF
ENDIF
;*****************************************

st_draw_options = (*Pst).opt

base_main  =WIDGET_BASE(/ROW, TITLE=title, MBAR=barBase, TLB_FRAME_ATTR=1, MAP=0,$
	BITMAP=image_bitmap)
widget_id = base_main

base_1     = WIDGET_BASE(base_main, /COLUMN, FRAME=0)
	base_rell = WIDGET_BASE(base_1, YSIZE=5)
	base_1a   = WIDGET_BASE(base_1, /COLUMN, FRAME=1)
	base_rell = WIDGET_BASE(base_1, YSIZE=5)
	base_1b   = WIDGET_BASE(base_1, /COLUMN, FRAME=1)
	base_rell = WIDGET_BASE(base_1, YSIZE=5)
	base_1c   = WIDGET_BASE(base_1, /COLUMN, FRAME=1)
	base_rell = WIDGET_BASE(base_1, YSIZE=5)
	base_1d   = WIDGET_BASE(base_1, /COLUMN, FRAME=0)
		base_1d1   = WIDGET_BASE(base_1, /COLUMN, /NONEXCLUSIVE)
base_2t    = WIDGET_BASE(base_main, /COLUMN, FRAME= 0)
	base_rell = WIDGET_BASE(base_2t, YSIZE=5)
	base_2    = WIDGET_BASE(base_2t, /COLUMN, FRAME= 1)
		base_rell = WIDGET_BASE(base_2, YSIZE=5)
		base_2a   = WIDGET_BASE(base_2, /COLUMN, FRAME=0)
		base_rell = WIDGET_BASE(base_2, YSIZE=5)
		base_2b   = WIDGET_BASE(base_2, /COLUMN, FRAME=0)
		base_2c   = WIDGET_BASE(base_2, /COLUMN, FRAME=0)

Filebar = WIDGET_BUTTON(barBase, VALUE='File', UVALUE='FILE', /MENU)
    Quitbar = WIDGET_BUTTON(Filebar, VALUE='Quit', UVALUE='BQUIT', /SEPARATOR)

ExportBar = WIDGET_BUTTON(barBase, VALUE='Export\Import', /MENU)

	exportImages  = WIDGET_BUTTON(Exportbar, VALUE='Export image',  UVALUE='EXPORT_IMAGE', UNAME='EXPORT_WINDOW', /MENU)
		exportbar_window  = WIDGET_BUTTON(ExportImages, VALUE='Export window',     UVALUE='EXPORT_IMAGE', UNAME='EXPORT_WINDOW')
		;exportbar_image   = WIDGET_BUTTON(ExportImages, VALUE='Export Image',      UVALUE='EXPORT_IMAGE', UNAME='EXPORT_IMAGE')
		exportbar_roi     = WIDGET_BUTTON(ExportImages, VALUE='Export ROI detail', UVALUE='EXPORT_IMAGE', UNAME='EXPORT_IMAGEROI')
	ExportParameters  = WIDGET_BUTTON(exportbar,     VALUE='Export ROI parameters',  /MENU)
		bttn_Exportmatrixformat = WIDGET_BUTTON(ExportParameters, VALUE='Matrix Format',  /SEPARATOR, UVALUE=['EXPORT_ASCII','MATRIX_FORMAT'], UNAME='EXPORT_ROI_RESULTS')
		bttn_Exportcolumnformat   = WIDGET_BUTTON(ExportParameters, VALUE='Single column format', UVALUE=['EXPORT_ASCII','COLUMN_FORMAT'], UNAME='EXPORT_ROI_RESULTS')
		bttn_Exportmultcolumnformat   = WIDGET_BUTTON(ExportParameters, VALUE='Multiple column format', 	UVALUE=['EXPORT_ASCII_ROI_MULTICOLUMN'], UNAME='EXPORT_ROI_MULTICOLUMN')
	bttn_ExportROIkinetics = WIDGET_BUTTON(exportbar, VALUE='Export ROI kinetic curve', UVALUE=['EXPORT_ASCII','COLUMN_FORMAT'], UNAME='EXPORT_KINETIC_CURVE', /SEPARATOR)
	bttn_ExportROIkinetics_Outliers = WIDGET_BUTTON(exportbar, VALUE='Export ROI kinetics (exclude Outliers)', UVALUE=['EXPORT_ASCII','COLUMN_FORMAT'], UNAME='EXPORT_KINETIC_CURVE_OUTLIERS')
	bttn_ExportRoi = WIDGET_BUTTON(exportbar, VALUE='Export ROI', /SEPARATOR,	UVALUE='EXPORT_ROI', UNAME='EXPORT_ROI_SAV')

Viewbar = WIDGET_BUTTON(barBase, VALUE='View', UVALUE='VIEW', /MENU)
    viewbar_ROIkinetics = WIDGET_BUTTON(viewbar, VALUE='View ROI kinetics', UVALUE='ROI_KINETICS', UNAME='VIEW_ROI_KINETICS', /SEPARATOR)
    viewbar_ROImodel = WIDGET_BUTTON(viewbar, VALUE='View ROI model', UVALUE='ROI_MODEL', UNAME='VIEW_ROI_MODEL')
    viewbar_ROIkinetics_NoOutliers = WIDGET_BUTTON(viewbar, VALUE='View ROI kinetics (exclude Outliers)', UVALUE='ROI_KINETICS', UNAME='VIEW_ROI_KINETICS_NOOUTLIERS', /SEPARATOR)
    ;viewbar_ROImodel_NoOutliers = WIDGET_BUTTON(viewbar, VALUE='View ROI model (exclude Outliers)', UVALUE='ROI_MODEL', UNAME='VIEW_ROI_MODEL_NOOUTLIERS')
    viewbar_deleteWindows = WIDGET_BUTTON(viewbar, VALUE='Close secondary windows', UVALUE='CLOSE_WINDOWS', /SEPARATOR)

Optionsbar  = WIDGET_BUTTON(barBase, VALUE='Options', UVALUE='OPTIONS', /MENU)

	barmenu_palettes = WIDGET_BUTTON(Optionsbar, VALUE='Color palette', /MENU, /SEPARATOR)
	    bar_palette1 = WIDGET_BUTTON(barmenu_palettes, VALUE='MRI', UVALUE='PALETTE1')
	    bar_palette3 = WIDGET_BUTTON(barmenu_palettes, VALUE='ROI', UVALUE='PALETTE3')
	bar_background_signal = WIDGET_BUTTON(Optionsbar, VALUE='Change background', 	UVALUE='BACKGROUND_IMAGE_SIGNAL')
	bar_colorbar_font     = WIDGET_BUTTON(Optionsbar, VALUE='Change colorbar font', 	UVALUE='COLORBAR_FONT', /MENU)
		bar_colorbar_mas  = WIDGET_BUTTON(bar_colorbar_font, VALUE='Increase size', 	UVALUE='COLORBAR_FONT', UNAME='INCREASE')
		bar_colorbar_menos= WIDGET_BUTTON(bar_colorbar_font, VALUE='Decrease size', 	UVALUE='COLORBAR_FONT', UNAME='DECREASE')
	bar_PlotsOpened_type = WIDGET_BUTTON(Optionsbar, VALUE='View plots', /MENU, /SEPARATOR)
		bttn_PlotsOpened_type1 =  WIDGET_BUTTON(bar_PlotsOpened_type, VALUE='Create new window', /CHECKED_MENU,	UVALUE='PLOTSOPENED_TYPE', UNAME='PLOTSOPENED_TYPE_1')
		bttn_PlotsOpened_type2 = WIDGET_BUTTON(bar_PlotsOpened_type, VALUE='Replace last window',  /CHECKED_MENU, UVALUE='PLOTSOPENED_TYPE', UNAME='PLOTSOPENED_TYPE_2')

Helpbar = WIDGET_BUTTON(barBase, VALUE='Help', UVALUE='BABOUT', /MENU)
	;bttn_usermanual = WIDGET_BUTTON(Helpbar, VALUE='User manual', UVALUE='BOPEN_USERMANUAL', SENSITIVE=0)
	bttn_about = WIDGET_BUTTON(Helpbar, VALUE='About DCE@urLAB...', UVALUE='BABOUT')

st_widgets_zoom     = fwBaseStruct_Zoom(PARENT_BASE=base_1a,MAIN_BASE=base_1a1, FRAME=0)
st_widgets_transp   = fwBaseStruct_Transparency(PARENT_BASE=base_1b, MAIN_BASE=base_1b1, FRAME=0)

st_widgets_parameter = fwBaseStruct_KineticParameter(PARENT_BASE=base_1c, MAIN_BASE=base_1c1, TITLE=(*(*Pst).modelinfo).str_model, VALUES=(*(*Pst).modelinfo).droplparams )

idx_paletteROI_ini   = 33l
idx_paletteImage_ini  = 0l

bttn_roilimits     = WIDGET_BUTTON(base_1d1, VALUE='Show ROI limits', UVALUE='ROI_LIMITS')

st_widgets_frames  = fwBaseStruct_Frames(PARENT_BASE=base_2a, MAIN_BASE=base_2a1, FRAME=0, $
	SET_VALUES=[(*Pst).inf.current_frame, (*Pst).inf.current_slice, (*Pst).inf.nframes, (*Pst).inf.nslices], ACTIVE=[1,0])

st_drawview  = fwBaseStruct_Window(PARENT_BASE=base_2b, TITLE='Dynamic image', MAIN_BASE=base_2b1, $
	ST_OPTIONS=st_draw_options, DRAW_ID=draw1, /ROI, UVALUE='DRAW', UNAME='WINDOW1', ST_OBJECTS=(*Pst).vi1)

strarr_text1 = STRARR(1) + '.'

text_info1  = WIDGET_TEXT(base_2c, VALUE=strarr_text1, /ALIGN_LEFT, XSIZE=30, YSIZE=1);, FONT=font1)

;---------------------------------------------------
st_widgets = { $
    wbase_main     : base_main,$
    wbases_resize1 : [base_1a1, base_1b1, base_1c1],$
    wbases_align   : [base_1a, base_2],$
    ;-----------------------------------------------
    wbttn_interp : st_widgets_zoom.bttn_interp,$
    ;-----------------------------------------------
    wlabel_n     : st_widgets_frames.label_n,$
    wlabel_z     : st_widgets_frames.label_z,$
    wslide_n     : st_widgets_frames.slide_n,$
    wslide_z     : st_widgets_frames.slide_z,$
    ;-----------------------------------------------
    wslide_transparency  : st_widgets_transp.slide_transparency,  $
    wlabel_transparency  : st_widgets_transp.label_transparency,  $
    wtext_palette_minval : st_widgets_transp.text_palette_minval, $
    wtext_palette_maxval : st_widgets_transp.text_palette_maxval, $
    Wdropl_roidrawmode   : st_widgets_transp.dropl_roidrawmode, $
    ;-----------------------------------------------
    wbttn_stddev         : st_widgets_parameter.bttn_stddev, $
    wbttn_value          : st_widgets_parameter.bttn_value,  $
    dropls_params        : st_widgets_parameter.dropls_params, $
    ;-----------------------------------------------
    bttn_roilimits       : bttn_roilimits,       $

    wbttn_Exportmatrixformat : bttn_Exportmatrixformat,$
    wbttn_Exportcolumnformat : bttn_Exportcolumnformat,$
    wbttn_Exportmultcolumnformat : bttn_Exportmultcolumnformat,$
    wbttn_ExportROIkinetics : bttn_ExportROIkinetics,$
    wbttn_ExportROIkinetics_outliers : bttn_ExportROIkinetics_outliers,$
    wbttn_PlotsOpened : [bttn_PlotsOpened_type1, bttn_PlotsOpened_type2], $
    wtext_info1          : text_info1,  $
    wdraw1               : draw1				 $
 }

;---------------------------------------------------

n_images_right = 3; (RCE, IAUC y TTM)

st_results = fStruct_Create(/RESULTSTRUCT, ST_IN=(*Pst).results, /COPY)
st_roi     = fStruct_Create(/ROISTRUCT,    ST_IN=(*Pst).roi, /COPY)

;---------------------------------------------------
sti={ $
	wbase_plot    : 0l,$
	wbase_parent  : (*Pst).wd.wbase_main,$
	data   : PTR_NEW(),    $        ; donde alojar los datos de la imagen dinamica T1
	mask   : PTR_NEW(),    $        ; mascara auxiliar
	parametric_AIF : PTR_NEW(),    $        ; parametric values of AIF (or Cp) are stored in a pointer, because can vary
	inf    : (*Pst).inf,                $
	opt    : (*Pst).opt,                $
	objplots : PTR_NEW(), $    ; Pointer to configuration of graphic objects (read from .sav file)

	wd  : st_widgets,                   $          ; Estructura con los identificadores de widgets
	vi1 : st_drawview,                  $          ; estructura con los objetos gráficos de dibujo de la imagen
	results   : TEMPORARY(st_results),  $
	roi       : TEMPORARY(st_roi),   $
	modelinfo : (*Pst).modelinfo,    $
	par       : (*Pst).par,          $
	flag_stddev : 0l,                $          ; Se pone a uno cuando se quiere trabajar con la varianza de la estimacion, no con la media
	flag_interp : (*Pst).flag_interp,    $
	flag_roilimits       : 0L, $
	opt_drawROIoutliers  : 0l, $  ; ROI outlayers (outside of limits of SCALE) are: 0) drawn 1) transparent 2) black
	opt_typeRCE          : 1l, $   ; tipo de calculo RCE (1, absoluto max/min, 2- relativo a first frame)
	opt_current_model    : (*Pst).opt_current_model,  $
	opt_current_param    : (*Pst).opt_current_param,  $
	opt_plotsOpened      : 2l, $, en la opcion por defecto, las ventanas machacan a las anteriores (2) en la otra, se van poniendo encima
	initial_zoom         : (*Pst).initial_zoom, $
	path_result          : (*Pst).path_result   $
}
;----------------------------------------------------------------------------------------------

Psti = PTR_NEW(Sti, /NO_COPY)

; Copy every data inside a pointer...
(*Psti).data      = PTR_NEW(*(*Pst).data)
(*Psti).mask      = PTR_NEW(*(*Pst).mask)
(*Psti).modelinfo = PTR_NEW(*(*Pst).modelinfo)

(*Psti).objplots  = PTR_NEW(*(*Pst).objplots)  ; copy content

ok  = SetWidget_maxsizes((*Psti).wd.wbases_resize1, /XSIZE)

ok  = fwdcemrirs_adjustScales(Psti)

LOADCT, 0
WIDGET_CONTROL, (*Psti).wd.wbase_main, /REALIZE, /MAP, SET_UVALUE=Psti

ok = SetWidget_RelativePosition((*Psti).wd.wbase_main, ID_PARENT=(*Psti).wbase_parent)

(*Psti).vi1.oimageROI -> setProperty, HIDE=0

ok = fwdcemri_Reshape_Windows(Psti)
ok = fwdcemri_Draw_Image(Psti, CONTROL_ZOOM=(*Psti).initial_zoom, NODRAW=1)
ok = fwdcemrirs_Redraw_im1(Psti)

WIDGET_CONTROL, (*Psti).wd.wdropl_roidrawmode, SET_DROPLIST_SELECT=(*Psti).opt_drawROIoutliers, /UPDATE
WIDGET_CONTROL, (*Psti).wd.dropls_params, SET_DROPLIST_SELECT=0, /UPDATE
WIDGET_CONTROL, (*Psti).wd.dropls_params, SEND_EVENT={ID:0L, TOP:0L, HANDLER:0L}

;-----------------------------------------------------------------------
XMANAGER, 'Interface_DCEMRIresult', (*Psti).wd.wbase_main, NO_BLOCK=1, GROUP_LEADER = (*Psti).wbase_parent, $
    EVENT_HANDLER='Interface_DCEMRIresult_event', CLEANUP='Interface_DCEMRIresult_Destroy'
;-----------------------------------------------------------------------

END ; Interface_DCEMRI

;****************************************************************************************************************************
;****************************************************************************************************************************
;<-->
;<++>

PRO Interface_InitialParameters_event, ev

WIDGET_CONTROL, ev.top, GET_UVALUE=id
WIDGET_CONTROL, ev.id,  GET_UVALUE=uval

opt_preview = 1

IF N_ELEMENTS(uval) GT 0 THEN BEGIN
CASE uval OF
'EXIT1':BEGIN

	WIDGET_CONTROL, id.wtext_Ah_ini,     GET_VALUE=str_Ah_ini
	WIDGET_CONTROL, id.wtext_keph_ini,   GET_VALUE=str_keph_ini
	WIDGET_CONTROL, id.wtext_kelh_ini,   GET_VALUE=str_kelh_ini

	WIDGET_CONTROL, id.wtext_ktrans_ini, GET_VALUE=str_ktrans_ini
	WIDGET_CONTROL, id.wtext_kep_ini,    GET_VALUE=str_kep_ini

	WIDGET_CONTROL, id.wtext_Sl_ini,      GET_VALUE=str_Sl_ini
	WIDGET_CONTROL, id.wtext_kepl_ini,    GET_VALUE=str_kepl_ini

	(*id.ptr_result).Ah_ini     = FLOAT(str_Ah_ini[0])
	(*id.ptr_result).keph_ini   = FLOAT(str_keph_ini[0]) ; En min-1
	(*id.ptr_result).kelh_ini   = FLOAT(str_kelh_ini[0]) ; En min-1

	(*id.ptr_result).ktrans_ini = FLOAT(str_ktrans_ini[0]) ; En min-1
	(*id.ptr_result).kep_ini    = FLOAT(str_kep_ini[0])    ; En min-1

	(*id.ptr_result).Sl_ini = FLOAT(str_Sl_ini[0]) ; En min-1
	(*id.ptr_result).kepl_ini    = FLOAT(str_kepl_ini[0])    ; En min-1
	(*id.ptr_result).ok = 1

	WIDGET_CONTROL, id.wbase_main, SET_UVALUE= id

	Interface_datain_destroy, ev.top

END
'EXIT2':BEGIN

	(*id.ptr_result).ok = 0
	WIDGET_CONTROL, id.wbase_main, SET_UVALUE= id

	Interface_datain_destroy, ev.top
END
'CHANGE VALUE: ': BEGIN
	; do nothing, all at exit1


END

ELSE:
ENDCASE
ENDIF

END

;*************************************************************************************
;*************************************************************************************

PRO Interface_InitialParameters_destroy, ev

	WIDGET_CONTROL, ev, GET_UVALUE=id
	WIDGET_CONTROL,  id.wbase_main,  /DESTROY

END
;*************************************************************************************
;*************************************************************************************

FUNCTION Interface_InitialParameters, GROUP_LEADER=group_leader, POSITION=position, $
    RESULT=result, INITIAL_RESULT=initial_result

title = 'Kinetics models: Initial parameters'

;tlb_frame_attr=1
tlb_frame_attr=9 ; dont close


IF N_ELEMENTS(group_leader) NE 0 THEN BEGIN
    base_floating = WIDGET_BASE(GROUP_LEADER=group_leader, FLOATING=1, UVALUE='INFO_FLOATING', $
       TITLE=title, TLB_FRAME_ATTR=tlb_frame_attr, MAP=opt_map, COLUMN=1)
ENDIF ELSE BEGIN
    base_floating = WIDGET_BASE(UVALUE='INFO_FLOATING', $
       TITLE=title, TLB_FRAME_ATTR=tlb_frame_attr, MAP=opt_map, COLUMN=1)
ENDELSE

base_0    = WIDGET_BASE(base_floating, COLUMN=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)

base_tab = WIDGET_TAB(base_0, /ALIGN_CENTER, LOCATION=0)
	base_t1 = WIDGET_BASE(base_tab, COLUMN=2, /ALIGN_LEFT, /BASE_ALIGN_LEFT, TITLE='Hoffmann model')
		base_t1a = WIDGET_BASE(base_t1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t1b = WIDGET_BASE(base_t1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t1c = WIDGET_BASE(base_t1, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
	base_t2 = WIDGET_BASE(base_tab, COLUMN=2, /ALIGN_LEFT, /BASE_ALIGN_LEFT, TITLE='Tofts model')
		base_t2a = WIDGET_BASE(base_t2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t2b = WIDGET_BASE(base_t2, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
	base_t3 = WIDGET_BASE(base_tab, COLUMN=2, /ALIGN_LEFT, /BASE_ALIGN_LEFT, TITLE='Larsson model')
		base_t3a = WIDGET_BASE(base_t3, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_t3b = WIDGET_BASE(base_t3, /ROW, /ALIGN_LEFT, /BASE_ALIGN_LEFT)


;IF N_ELEMENTS(initial_result) EQ 0 THEN RETURN, -1
IF N_ELEMENTS(initial_result) EQ 0 THEN BEGIN

	initial_result = {ktrans_ini: 0d,$
	   kep_ini :  0d,  $
	   Ah_ini  :  0d,  $
	   keph_ini:  0d,$
	   kelh_ini:  0d,$
	   sl_ini  :  0d,$
	   kepl_ini:  0d,$
	   ok : 0b $
	 }
ENDIF

str_ktrans_ini  = STRTRIM(STRING(initial_result.ktrans_ini, FORMAT='(F10.5)'),2) ; EN min-1
str_kep_ini     = STRTRIM(STRING(initial_result.kep_ini,    FORMAT='(F10.5)'),2) ; en min-1

str_Ah_ini    = STRTRIM(STRING(initial_result.Ah_ini, FORMAT='(F10.5)'),2) ; sin unidades
str_keph_ini  = STRTRIM(STRING(initial_result.keph_ini,    FORMAT='(F10.5)'),2) ; en min-1
str_kelh_ini  = STRTRIM(STRING(initial_result.kelh_ini,    FORMAT='(F10.5)'),2) ; en min-1

str_Sl_ini    = STRTRIM(STRING(initial_result.Sl_ini, FORMAT='(F10.5)'),2) ; sin unidades
str_kepl_ini  = STRTRIM(STRING(initial_result.kepl_ini,    FORMAT='(F10.5)'),2) ; en min-1



text_Ah_ini  = WIDGET_TEXT(base_t1a, VALUE=str_Ah_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_Ah_ini = WIDGET_LABEL(base_t1a, VALUE=' A (n.d)')

text_keph_ini    = WIDGET_TEXT(base_t1b, VALUE=str_keph_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1,  /EDITABLE, UVALUE=uvalue)
label_keph_ini   = WIDGET_LABEL(base_t1b, VALUE=' Kep(h) (min-1)')

text_kelh_ini    = WIDGET_TEXT(base_t1c, VALUE=str_kelh_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1,  /EDITABLE, UVALUE=uvalue)
label_kelh_ini   = WIDGET_LABEL(base_t1c, VALUE=' Kel(h) (min-1)')

text_ktrans_ini = WIDGET_TEXT(base_t2a, VALUE=str_ktrans_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_ktrans_ini = WIDGET_LABEL(base_t2a, VALUE=' KTRANS  (min-1)')

text_kep_ini    = WIDGET_TEXT(base_t2b, VALUE=str_kep_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1,  /EDITABLE, UVALUE=uvalue)
label_kep_ini   = WIDGET_LABEL(base_t2b, VALUE=' Kep  (min-1)')

text_Sl_ini = WIDGET_TEXT(base_t3a, VALUE=str_Sl_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1, /EDITABLE, UVALUE=uvalue)
label_Sl_ini = WIDGET_LABEL(base_t3a, VALUE=' S(l)  (n,d)')

text_kepl_ini    = WIDGET_TEXT(base_t3b, VALUE=str_kep_ini, XSIZE=8, $
	/ALIGN_LEFT, SENSITIVE=1,  /EDITABLE, UVALUE=uvalue)
label_kepl_ini   = WIDGET_LABEL(base_t3b, VALUE=' Kep(l)  (min-1)')


base_1 = WIDGET_BASE(base_floating, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
	bttn_close  = WIDGET_BUTTON(base_1, VALUE=' OK     ', UVALUE='EXIT1')
	bttn_cancel = WIDGET_BUTTON(base_1, VALUE=' Cancel ', UVALUE='EXIT2')

st = { $
	wbase_main     :   base_floating,  $
	wbttn_close    :   bttn_close,     $
	wtext_Ah_ini   :   text_Ah_ini,    $
	wtext_keph_ini :   text_keph_ini,  $
	wtext_kelh_ini :   text_kelh_ini,  $
	wtext_ktrans_ini : text_ktrans_ini,$
	wtext_kep_ini :    text_kep_ini,   $
	wtext_kepl_ini :   text_kepl_ini,   $
	wtext_Sl_ini :    text_Sl_ini,   $
	ptr_result  : PTR_NEW(initial_result, /ALLOCATE_HEAP),$
	last : 0l }


IF N_ELEMENTS(position) EQ 1 AND N_ELEMENTS(group_leader) NE 0 THEN BEGIN
	pos = GetWidget_RelativePosition(base_floating, ID_PARENT=group_leader, POSITION=position)
	WIDGET_CONTROL, base_floating, XOFFSET=pos[0], YOFFSET=pos[1]
ENDIF ELSE BEGIN
	 pos = GetWidget_PositionInScreen(base_floating,CENTER=1)
ENDELSE

;-----------------------------------------------------------------

WIDGET_CONTROL, base_t3, MAP=0


WIDGET_CONTROL, base_floating, SET_UVALUE=st
WIDGET_CONTROL, base_floating, /REALIZE

XMANAGER,  'Interface_InitialParameters',  base_floating, NO_BLOCK=1, MODAL=1, $
    EVENT_HANDLER='Interface_InitialParameters_event', CLEANUP='Interface_InitialParameters_Destroy'

PRINT, 'OK'

result = *st.ptr_result

PTR_FREE, st.ptr_result
; Interfaz modal de entrada de datos

RETURN, result.ok

END

;<-->
;<+>
FUNCTION Interface_LabelInfo, GROUP_LEADER=group_leader, XSIZE=xsize, YSIZE=ysize, POSITION=position,$
    TITLE=title, LABEL=label, NO_MAP=no_map, NO_CLOSE=no_close, NOTITLE=notitle, MODAL=modal, $
    IMAG = imag, FONT_SIZE=font_size, FONT_TYPE=font_type, FONT_EFFECT=font_effect, IMAGE_LOGO=image_logo,$
    HORIZONTAL=horironzal, OKCANCEL=okcancel

IF KEYWORD_SET(horizontal) THEN BEGIN
	opt_horizontal=1
	opt_vertical = 0
	xoffset_label = 20
	yoffset_label = 10
ENDIF ELSE BEGIN
	opt_horizontal=0
	opt_vertical = 1
	xoffset_label = 10
	yoffset_label = 20
ENDELSE

; la fuente por defecto es 12, arial

n_labels = N_ELEMENTS(label)
IF n_labels EQ 0 THEN RETURN, -1
IF SIZE(label, /TYPE) NE 7 THEN RETURN, -1

CASE N_ELEMENTS(font_size) OF
0 : arr_font_size = INTARR(n_labels)+14
1 : arr_font_size = REPLICATE(font_size, n_labels)
n_labels: arr_font_size = FIX(font_size)
ELSE: RETURN,-1
ENDCASE

CASE N_ELEMENTS(font_type) OF
0 : arr_font_type = STRARR(n_labels)+'arial'
1 : arr_font_type = REPLICATE(font_type, n_labels)
n_labels: arr_font_type = STRING(font_type)
ELSE: RETURN,-1
ENDCASE

CASE N_ELEMENTS(font_effect) OF
0 : arr_font_effect = STRARR(n_labels) + 'light'
1 : arr_font_effect = REPLICATE(font_effect, n_labels)
n_labels: arr_font_effect = STRING(font_effect)
ELSE: RETURN,-1
ENDCASE

; Crea una base sencilla para poner información (about...) y opcionalmente una imagen

; POSICION: Array de dos elementos indicando los pixeles en X e Y de la esquina superior derecha
; por defecto, en el centro.

IF N_ELEMENTS(title) EQ 0 THEN tit = '' ELSE tit=STRING(title[0])
IF KEYWORD_SET(no_map) THEN opt_map = 0 ELSE opt_map = 1
IF KEYWORD_SET(no_close) THEN tlb_frame_attr=11+16L ELSE tlb_frame_attr=1+16L

IF N_ELEMENTS(group_leader) NE 0 THEN BEGIN
    base_floating = WIDGET_BASE(GROUP_LEADER=group_leader, FLOATING=1, UVALUE='INFO_FLOATING', $
       TITLE=tit, TLB_FRAME_ATTR=tlb_frame_attr, MAP=opt_map, COLUMN=1, MODAL=modal)
ENDIF ELSE BEGIN
    base_floating = WIDGET_BASE(UVALUE='INFO_FLOATING', $
       TITLE=tit, TLB_FRAME_ATTR=tlb_frame_attr, MAP=opt_map, COLUMN=1)
ENDELSE

base_2   = WIDGET_BASE(base_floating,  ROW=opt_horizontal, COLUMN=opt_vertical, $
		/ALIGN_CENTER, /BASE_ALIGN_CENTER)
	base_2a = WIDGET_BASE(base_2,  ROW=1, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
	base_2b = WIDGET_BASE(base_2,  ROW=1, /ALIGN_CENTER, /BASE_ALIGN_CENTER)

opt_im = 0l
IF N_ELEMENTS(image_logo) GT 3 THEN BEGIN ; ademas añade un logo
	image_size = SIZE(image_logo, /DIMENSIONS)
	IF image_size[0] EQ 3 AND N_ELEMENTS(image_size) EQ 3 THEN BEGIN
		opt_im = 1
	ENDIF
ENDIF

IF opt_im EQ 1 THEN BEGIN
	wdraw_logo = WIDGET_DRAW(base_2a,  COLOR_MODEL=0, FRAME=0, XSIZE=image_size[1], YSIZE=image_size[2])
ENDIF

base_2a    = WIDGET_BASE(base_2B, XSIZE=xoffset_label, YSIZE=yoffset_label)
base_info  = WIDGET_BASE(base_2B, COLUMN=1)
base_2b    = WIDGET_BASE(base_2B, XSIZE=xoffset_label, YSIZE=yoffset_label, /ALIGN_RIGHT)

FOR i=0l, n_labels-1 DO BEGIN
	opt_font = STRTRIM(arr_font_type[i],2) + '*' + arr_font_effect[i] + '*' + STRTRIM(arr_font_size[i],2)
	label_info = WIDGET_LABEL(base_info, VALUE=label[i], FONT=opt_font, /ALIGN_LEFT)
ENDFOR

IF KEYWORD_SET(okcancel) THEN BEGIN
	base_2c    = WIDGET_BASE(base_2, /ROW, /ALIGN_RIGHT)
	button_ok     = WIDGET_BUTTON(base_2c, VALUE='OK',UVALUE='OK',    EVENT_FUNC='Interface_LabelInfo_ok')
	base_rell     = WIDGET_BASE(base_2c, XSIZE = 10)
	button_cancel = WIDGET_BUTTON(base_2c, VALUE='Cancel',UVALUE='CANCEL',EVENT_FUNC='Interface_LabelInfo_cancel')
	base_rell     = WIDGET_BASE(base_2C, XSIZE = 10)
ENDIF

IF N_ELEMENTS(position) EQ 2 THEN pos = position
IF N_ELEMENTS(position) EQ 0 AND N_ELEMENTS(group_leader) NE 0 THEN BEGIN
	pos = GetWidget_RelativePosition(base_floating, ID_PARENT=group_leader, POSITION=4)
ENDIF

IF N_ELEMENTS(pos) EQ 0 THEN BEGIN
	ok = GetWidget_PositionInScreen(base_floating, CENTER=1)
ENDIF ELSE BEGIN
	WIDGET_CONTROL, base_floating, XOFFSET=pos[0], YOFFSET=pos[1]
ENDELSE

WIDGET_CONTROL, base_floating, /REALIZE, MAP=1

IF opt_im EQ 1 THEN BEGIN
	WIDGET_CONTROL, wdraw_logo,GET_VALUE=wdraw_id
	WSET, wdraw_id
	TV, image_logo, TRUE=1
	WSET, -1
ENDIF

IF KEYWORD_SET(okcancel) THEN BEGIN
	ptr_result = PTR_NEW('CANCEL')
	WIDGET_CONTROL, base_floating, SET_UVALUE=ptr_result
	XMANAGER,  'Interface_labelinfo',  base_floating, NO_BLOCK=1, MODAL=1
	RETURN, *ptr_result

ENDIF ELSE BEGIN
	RETURN, base_floating
ENDELSE

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<++>

FUNCTION Get_buttonsIcons, PLAY=play, PAUSE=pause, STOPI=stopi, RIGHT=right, LEFT=left, ROTATE_CLOCKWISE=rotate_clockwise, ROTATE_ANTICLOCKWISE=rotate_anticlockwise, $
	FAST_REWIND_RIGHT=fast_rewind_right, FAST_REWIND_LEFT=fast_rewind_left

IF KEYWORD_SET(play) THEN BEGIN
var = BYTE([ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 1, 0], [0, 7, 0], [0, 31, 0], [0, 127, 0],  $
[0, 255, 0], [0, 255, 3], [0, 255, 1], [0, 127, 0], [0, 31, 0], [0, 15, 0], [0, 3, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],  $
[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] )
ENDIF
IF KEYWORD_SET(pause) THEN BEGIN
var = BYTE([ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 199, 1], [0, 199, 1], [0, 199, 1], [0, 199, 1],  $
[0, 199, 1], [0, 199, 1], [0, 199, 1], [0, 199, 1], [0, 199, 1], [0, 199, 1], [0, 199, 1], [0, 0, 0], [0, 0, 0], [0, 0, 0],  $
[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] )
ENDIF
IF KEYWORD_SET(stopi) THEN BEGIN
var = BYTE([ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [128, 255, 3], [128, 255, 3], [128, 255, 3], [128, 255, 3],  $
[128, 255, 3], [128, 255, 3], [128, 255, 3], [128, 255, 3], [128, 255, 3], [128, 255, 3], [128, 255, 3], [0, 0, 0], [0, 0, 0], [0, 0, 0],  $
[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] )
ENDIF
IF KEYWORD_SET(right) THEN BEGIN
var = BYTE([ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 24, 0], [0, 120, 0], [0, 248, 1],  $
[248, 255, 7], [248, 255, 31], [248, 255, 31], [248, 255, 7], [0, 248, 1], [0, 120, 0], [0, 24, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],  $
[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] )
ENDIF
IF KEYWORD_SET(left) THEN BEGIN
var = BYTE([ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 24, 0], [0, 30, 0], [128, 31, 0],  $
[224, 255, 31], [248, 255, 31], [248, 255, 31], [224, 255, 31], [128, 31, 0], [0, 30, 0], [0, 24, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],  $
[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] )
ENDIF
IF KEYWORD_SET(rotate_clockwise) THEN BEGIN
var = BYTE([ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 24, 0], [0, 56, 0], [0, 127, 0], [128, 255, 0], [192, 127, 0], [224, 57, 0],  $
[224, 8, 0], [112, 0, 0], [112, 0, 0], [112, 0, 7], [112, 0, 7], [224, 128, 3], [224, 193, 3], [192, 255, 1], [192, 255, 0], [0, 127, 0],  $
[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] )
ENDIF
IF KEYWORD_SET(rotate_anticlockwise) THEN BEGIN
var = BYTE([ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 16, 0], [0, 24, 0], [0, 252, 0], [0, 254, 3], [0, 254, 7], [0, 152, 7],  $
[0, 16, 15], [0, 0, 14], [0, 0, 14], [192, 0, 14], [224, 0, 14], [192, 1, 15], [192, 3, 7], [128, 255, 7], [0, 255, 3], [0, 254, 0],  $
[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] )
ENDIF
IF KEYWORD_SET(fast_rewind_right) THEN BEGIN
	var = BYTE([ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 16, 0], [0, 50, 0], [0, 114, 0], [0, 246, 0],  $
	[0, 254, 1], [0, 254, 3], [0, 254, 3], [0, 246, 1], [0, 114, 0], [0, 50, 0], [0, 16, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],  $
	[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] )
ENDIF
IF KEYWORD_SET(fast_rewind_left) THEN BEGIN
	var = BYTE([ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 136, 0], [0, 204, 0], [0, 206, 0],  $
	[128, 239, 0], [192, 239, 0], [128, 239, 0], [0, 207, 0], [0, 206, 0], [0, 140, 0], [0, 8, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],  $
	[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] )
END

RETURN, var

END

;**************************************************************************************************
;**************************************************************************************************

FUNCTION fwPlayer_DisplayImage, Pst

	WIDGET_CONTROL,(*Pst).wd.slide_n,  GET_VALUE=n_n

	tam = SIZE((*Pst).data_bytscl, /DIMENSIONS)
    img_bytscl = REFORM(((*Pst).data_bytscl)[*,*,n_n])

    (*Pst).vi1.omodel  -> setProperty, HIDE=0
    (*Pst).vi1.oimage  -> setProperty, INTERPOLATE=(*Pst).flag_interp

	;min_val =  MIN(img,    MAX=max_val)
    ;img_scl =  BYTSCL(img, MIN=min_val, MAX=max_val)

    xMax = tam[0] & yMax = tam[1]

    (*Pst).vi1.oimage -> setProperty, DATA = img_bytscl, DIMENSIONS=[xMax, yMax], GREYSCALE=0, HIDE=0, SUB_RECT=[0,0,xMax,yMax]
    (*Pst).vi1.oview  -> setproperty, VIEWPLANE_RECT=[0,0,xMax, yMax]


	WIDGET_CONTROL, (*Pst).wd.Draw, GET_VALUE = window_1
    window_1 -> Draw, (*Pst).vi1.oView

	RETURN, 1

END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_Player_destroy ,event


END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_Player_event, event

WIDGET_CONTROL,event.top, GET_UVALUE = Pst
WIDGET_CONTROL,event.id,  GET_UVALUE = uvalue; valor que viene de cada evento

IF N_ELEMENTS(uvalue) EQ 0 THEN RETURN

;----------------------------------------------------------------------------------------
CASE uvalue OF

'PLAYER' : BEGIN

	n_frames = (*Pst).nframes;

	uname = WIDGET_INFO(event.id, /UNAME)
	CASE uname OF
	'START': BEGIN
		(*(*Pst).player).stop  = 0b
		(*(*Pst).player).pause = 0b
		WIDGET_CONTROL, (*Pst).wd.slide_n, GET_VALUE=frame, SENSITIVE=0
		(*Pst).current_frame  = frame
		(*(*Pst).player).iter = frame
		;WIDGET_CONTROL,(*(*Pst).player).wtime, TIMER=(*(*Pst).player).frame_defaulttime
		WIDGET_CONTROL,(*(*Pst).player).wtime, TIMER=(*(*Pst).player).time
		WIDGET_CONTROL, (*Pst).wd.label_info, SET_VALUE= STRTRIM(STRING(  (*(*Pst).player).time, FORMAT='(F10.5)'),2) + ' s/frame  '
	END
	'STOP': BEGIN
		(*(*Pst).player).stop = 1b
		WIDGET_CONTROL, (*Pst).wd.slide_n, SET_VALUE=0, SENSITIVE=1
		WIDGET_CONTROL, (*Pst).wd.label_n, SET_VALUE='0'
		(*Pst).current_frame = 0
		image = ((*Pst).data_bytscl)[*,*,0] ; already scaled

		(*Pst).vi1.oimage->setproperty, DATA=image
	   	(*Pst).vi1.oWindow->Draw,(*Pst).vi1.oview
		END
	'PAUSE': BEGIN
		(*(*Pst).player).pause =1b
		WIDGET_CONTROL, (*Pst).wd.slide_n, SENSITIVE=1
	END
	'FAST_REWIND_RIGHT': BEGIN ;accelerates the cine
		(*(*Pst).player).time = (((*(*Pst).player).time)/((*(*Pst).player).mult_inc) > ((*(*Pst).player).frame_mintime))
		WIDGET_CONTROL,(*(*Pst).player).wtime, TIMER=(*(*Pst).player).time
		WIDGET_CONTROL, (*Pst).wd.label_info, SET_VALUE= STRTRIM(STRING(  (*(*Pst).player).time, FORMAT='(F10.5)'),2) + ' s/frame  '

	END
	'FAST_REWIND_LEFT': BEGIN ;accelerates the cine
		(*(*Pst).player).time = (((*(*Pst).player).time)*((*(*Pst).player).mult_inc) < ((*(*Pst).player).frame_maxtime))
		WIDGET_CONTROL,(*(*Pst).player).wtime, TIMER=(*(*Pst).player).time
		WIDGET_CONTROL, (*Pst).wd.label_info, SET_VALUE= STRTRIM(STRING(  (*(*Pst).player).time, FORMAT='(F10.5)'),2) + ' s/frame  '
	END
	'TIMER': BEGIN
		;-----------------------------------------------------------------------
		i = (*(*Pst).player).iter MOD n_frames
		image = ((*Pst).data_bytscl)[*,*,i]
		(*Pst).vi1.oimage->setproperty, DATA=image
		(*Pst).vi1.oWindow->Draw,(*Pst).vi1.oview
		;-----------------------------------------------------------------------
		WIDGET_CONTROL, (*Pst).wd.slide_n, SET_VALUE=i
		WIDGET_CONTROL, (*Pst).wd.label_n, SET_VALUE=STRTRIM(i,2)

		;more to end the loop
		IF (*(*Pst).player).stop THEN BEGIN
			WIDGET_CONTROL, (*Pst).wd.slide_n, SET_VALUE=0
			(*Pst).current_frame = 0
			RETURN
		ENDIF
		IF (*(*Pst).player).pause THEN BEGIN
			(*Pst).vi1.oWindow->Draw,(*Pst).vi1.oview
			WIDGET_CONTROL, (*Pst).wd.slide_n, GET_VALUE=frame
			(*Pst).current_frame = frame
			RETURN
		ENDIF
		(*(*Pst).player).iter++
		WIDGET_CONTROL, (*(*Pst).player).wtime, TIMER=(*(*Pst).player).time
	END
	ENDCASE
END
;----------------------------------------------------------------------------------------
'S_CHANGE' : BEGIN
	WIDGET_CONTROL, (*Pst).wd.slide_n, GET_VALUE=frame
	image = ((*Pst).data_bytscl)[*,*,frame] ; already scaled
	(*Pst).vi1.oimage->setproperty, DATA=image
	(*Pst).vi1.oWindow->Draw,(*Pst).vi1.oview
END
;----------------------------------------------------------------------------------------
ELSE:
ENDCASE

END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_Player, DATA=data, PALETTE=opalette, SIZE_VIEW=size_view, SIZE_SCROLL=size_scroll

base_main    = WIDGET_BASE(/COLUMN, FRAME=0, /ALIGN_CENTER)

IF N_ELEMENTS(title) EQ 0 THEN title=''

IF N_ELEMENTS(size_view)   EQ 0 THEN size_view = [256,256l]*2
IF N_ELEMENTS(size_scroll) EQ 0 THEN size_scroll = [256,256l]*2


size_image = (SIZE(data, /DIMENSIONS))[0:2]   ; 3d image (only dynamic in one slice)
n_frames   = size_image[2]

current_frame = 0l
value_label_n = STRTRIM(STRING(current_frame),2)
data_bytscl   = BYTARR(size_image)

min_val =  MIN(data,   MAX=max_val) ; general scaling

FOR i=0l, n_frames-1 DO BEGIN
	img = data[*,*,i]
	data_bytscl[*,*,i] = BYTSCL(img, MIN=min_val, MAX=max_val)
ENDFOR

viewsize    = DBLARR(2)
scrollsize = DBLARR(2)
scrollsize[0] = size_scroll[0]
scrollsize[1] = (1d)*size_scroll[0]*size_image[1]/size_image[0]  ; new to avoid the palette

relative_change = scrollsize[1]/size_scroll[1];    different to 1 due to palette

viewsize[0] = size_view[0]
viewsize[1] = size_view[1]*relative_change



base_1      = WIDGET_BASE(base_main, /COLUMN, FRAME=0, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
	base_1a  = WIDGET_BASE(base_1, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
		base_1a1  = WIDGET_BASE(base_1a, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
	base_2a  = WIDGET_BASE(base_1, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
		base_player  = WIDGET_BASE(base_2a, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER, /FRAME)
		base_text   =  WIDGET_BASE(base_2a, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)

base_draw   = WIDGET_BASE(base_main, /COLUMN, FRAME=0, /ALIGN_CENTER, /BASE_ALIGN_CENTER)
draw        = WIDGET_DRAW(Base_draw, RETAIN=2, COLOR_MODEL=0, RENDERER=1,GRAPHICS_LEVEL=2, $
					SCR_XSIZE=viewsize[0], SCR_YSIZE=viewsize[1], $
					X_SCROLL_SIZE=scrollsize[0], Y_SCROLL_SIZE=scrollsize[1], SCROLL=1, $
					MOTION_EVENTS=0, BUTTON_EVENTS=0, TRACKING_EVENTS=0, VIEWPORT_EVENTS=1, EXPOSE_EVENTS=0,$
					UVALUE=uvalue, UNAME=uname)

slide_n = WIDGET_SLIDER(base_1a1, MAXIMUM=n_frames, MINIMUM=0,  SCR_XSIZE=size_slider_1, $
    VALUE=current_frame, UVALUE='S_CHANGE', UNAME = 'SLICE_N',  SENSITIVE=1, TITLE='Dynamic frame', /SUPPRESS_VALUE)
wlabel_n = WIDGET_LABEL(base_1a1, VALUE=value_label_n, SCR_XSIZE=25, SCR_YSIZE=25)
main_base = TEMPORARY(base_1)

size_rell = 8

button_size = [24,24l]
wbaseCine  = WIDGET_BASE(base_player,    UVALUE='PLAYER', UNAME='TIMER')
wrell       = WIDGET_BASE(base_player, XSIZE=size_rell)
wplayCine   = WIDGET_BUTTON(base_player,  EVENT_PRO='Interface_Player_event',$
	VALUE= get_buttonsIcons(/PLAY),   UVALUE='PLAYER', UNAME= 'START',/ALIGN_CENTER,XSIZE=button_size[0],YSIZE=button_size[1])
wrell       = WIDGET_BASE(base_player, XSIZE=size_rell)
wpauseCine  = WIDGET_BUTTON(base_player,  EVENT_PRO='Interface_Player_event',$
	VALUE= get_buttonsIcons(/PAUSE),  UVALUE='PLAYER',  UNAME='PAUSE',/ALIGN_CENTER,XSIZE=button_size[0],YSIZE=button_size[1])
wrell       = WIDGET_BASE(base_player, XSIZE=size_rell)
wstopCine   = WIDGET_BUTTON(base_player,  EVENT_PRO='Interface_Player_event',$
	VALUE= get_buttonsIcons(/STOPI),  UVALUE='PLAYER',  UNAME='STOP',/ALIGN_CENTER,XSIZE=button_size[0],YSIZE=button_size[1])
wrell       = WIDGET_BASE(base_player, XSIZE=size_rell)
wrewindleftCine   = WIDGET_BUTTON(base_player,  EVENT_PRO='Interface_Player_event',$
	VALUE= get_buttonsIcons(/FAST_REWIND_LEFT),  UVALUE='PLAYER',  UNAME='FAST_REWIND_LEFT',/ALIGN_CENTER,XSIZE=button_size[0],YSIZE=button_size[1])
wrell       = WIDGET_BASE(base_player, XSIZE=size_rell)
wrewindrightCine   = WIDGET_BUTTON(base_player,  EVENT_PRO='Interface_Player_event',$
	VALUE= get_buttonsIcons(/FAST_REWIND_RIGHT),  UVALUE='PLAYER',  UNAME='FAST_REWIND_RIGHT',/ALIGN_CENTER,XSIZE=button_size[0],YSIZE=button_size[1])
wrell       = WIDGET_BASE(base_player, XSIZE=size_rell)

wlabel_info  = WIDGET_LABEL(base_text, UVALUE='TEXT_INFO', VALUE='', XSIZE=120)


WIDGET_CONTROL, draw, GET_VALUE  = ds_window

IF N_ELEMENTS(opalette) NE 0 THEN BEGIN
	opalette->GetProperty, RED=r1,GREEN=g1,BLUE=b1
ENDIF ELSE BEGIN
	LOADCT,0
	TVLCT, r1,g1,b1, /GET
	opalette = OBJ_NEW('IDLgrPalette',r1,g1,b1)
ENDELSE

sto = {$
	owindow    : ds_window,              $
	oview      : OBJ_NEW('IDLgrView'),   $
	omodel     : OBJ_NEW('IDLgrModel'),  $
	oimage     : OBJ_NEW('IDLgrImage', PALETTE=opalette, GREYSCALE=0) $
}

sto.omodel ->add, sto.oimage
sto.oview  ->add, sto.omodel

st_widget = { $
	base_main: base_main,$
	slide_n : slide_n,$
	playCine : wplaycine,  $
	stopCine : wstopCine,  $
	pauseCine: WpauseCine, $
	draw     : draw, $
	label_info : wlabel_info, $
	label_n    : wlabel_n     $
}

st_player = {wtime:wbaseCine,$
		time   :0.128D, $
		stop   :0b,  $
		pause  :0b,  $
		frame_defaulttime : 0.128D, $
		frame_maxtime     : 1.024D, $ ; two seconds
		frame_mintime     : 0.016D, $
		mult_inc          : 2d, $
		iter:0l}
player = PTR_NEW(st_player,/NO_COPY)

Sti = {  $
    wd :st_widget,    $
	vi1:sto,          $
	player : player,  $
	flag_interp : 0l, $
	data_bytscl:data_bytscl ,$
	nframes:n_frames, $
	current_frame:current_frame $
}

Psti = PTR_NEW(Sti, /NO_COPY)

WIDGET_CONTROL, (*Psti).wd.base_main, /REALIZE, /MAP, SET_UVALUE=Psti

WIDGET_CONTROL, (*Psti).wd.label_info, SET_VALUE= STRTRIM(STRING((*(*Psti).player).time, FORMAT='(F10.5)'),2) + ' s/frame  '
WIDGET_CONTROL, (*Psti).wd.draw,  SCR_XSIZE=viewsize[0], SCR_YSIZE=viewsize[1], DRAW_XSIZE=scrollsize[0],  DRAW_YSIZE=scrollsize[1]


(*Psti).vi1.oview  -> setProperty,  DIMENSIONS=scrollsize, UNITS=0

WIDGET_CONTROL, (*Psti).wd.slide_n, SCR_XSIZE=viewsize[0]*0.60
WIDGET_CONTROL, (*Psti).wd.draw, GET_VALUE  = ds_window
(*Psti).vi1.owindow = ds_window

ok = fwPlayer_DisplayImage(Psti)

;-----------------------------------------------------------------------
XMANAGER, 'Interface_Player', (*Psti).wd.base_main, NO_BLOCK=1,$
    EVENT_HANDLER='Interface_Player_Event', CLEANUP='Interface_Player_Destroy'; $
    ;GROUP_LEADER = (*Psti).wbase_parent, $
;-----------------------------------------------------------------------

END

;**************************************************************************************************
;**************************************************************************************************

PRO Test_Interface_Player


file_sav = 'F:\Experiments_DCE-MRI\Tests\dce-mri_data.SAV'

RESTORE, file_sav

HELP, data

datadyn = REFORM(data[*,*,0,*])

Interface_Player, DATA=datadyn, SIZE_VIEW=size_view, SIZE_SCROLL=size_scroll

END

;<-->
;<++>

FUNCTION fwidget_ReadRaw_Event, id

	WIDGET_CONTROL, id.wtext_xsize, GET_VALUE=str_xsize
	WIDGET_CONTROL, id.wtext_ysize, GET_VALUE=str_ysize
	WIDGET_CONTROL, id.wtext_nslices, GET_VALUE=str_nslices
	WIDGET_CONTROL, id.wtext_nframes, GET_VALUE=str_nframes
	WIDGET_CONTROL, id.wtext_offset, GET_VALUE=str_offset

	WIDGET_CONTROL, id.wtext_injection, GET_VALUE=str_frame_injection
	WIDGET_CONTROL, id.wtext_frame_period, GET_VALUE=str_frame_period

	xsize = Get_IntegerFromString(str_xsize[0], NEW_STRING=str_xsize2, ERROR_ID=error_id, /POSITIVE, MAXVAL=1024)
	IF error_id EQ 0 THEN BEGIN
		(*id.ptr_st).xsize = xsize
	ENDIF
	WIDGET_CONTROL, id.wtext_xsize, SET_VALUE=str_xsize2

	ysize = Get_IntegerFromString(str_ysize[0], NEW_STRING=str_ysize2, ERROR_ID=error_id, /POSITIVE, MAXVAL=1024)
	IF error_id EQ 0 THEN BEGIN
		(*id.ptr_st).ysize = Ysize
	ENDIF
	WIDGET_CONTROL, id.wtext_ysize, SET_VALUE=str_ysize2

	nslices = Get_IntegerFromString(str_nslices[0], NEW_STRING=str_nslices2, ERROR_ID=error_id, /POSITIVE, MAXVAL=512)
	IF error_id EQ 0 THEN BEGIN
		(*id.ptr_st).nslices = nslices
	ENDIF
	WIDGET_CONTROL, id.wtext_nslices, SET_VALUE=str_nslices2

	nframes = Get_IntegerFromString(str_nframes[0], NEW_STRING=str_nframes2, ERROR_ID=error_id, /POSITIVE, MAXVAL=5000)
	IF error_id EQ 0 THEN BEGIN
		(*id.ptr_st).nframes = nframes
	ENDIF
	WIDGET_CONTROL, id.wtext_nframes, SET_VALUE=str_nframes2

	offset = Get_IntegerFromString(str_offset[0], NEW_STRING=str_offset2, ERROR_ID=error_id, /POSORZERO)
	IF error_id EQ 0 THEN BEGIN
		(*id.ptr_st).offset = offset
	ENDIF
	WIDGET_CONTROL, id.wtext_offset, SET_VALUE=str_offset2

	frame_injection = Get_IntegerFromString(str_frame_injection[0], NEW_STRING=str_frame_injection2, $
		ERROR_ID=error_id, /POSITIVE)
	IF error_id EQ 0 THEN BEGIN
		(*id.ptr_st).frame_injection =(frame_injection-1) > 0
	ENDIF
	WIDGET_CONTROL, id.wtext_injection, SET_VALUE=str_frame_injection2

	frame_period = Get_FloatingFromString(str_frame_period[0], NEW_STRING=str_frame_period2, $
		ERROR_ID=error_id, /POSITIVE, MAX_DECIMALS=2)
	IF error_id EQ 0 THEN BEGIN
		(*id.ptr_st).frame_period =frame_period
	ENDIF
	WIDGET_CONTROL, id.wtext_frame_period, SET_VALUE=str_frame_period2



	FOR i=0l, 3L DO BEGIN
		(*id.ptr_st).invarr[i] = WIDGET_INFO(id.wbttn_inv[i], /BUTTON_SET)
	ENDFOR
	(*id.ptr_st).littleendian = WIDGET_INFO(id.wbttn_littlendian, /BUTTON_SET)

	type_select      = WIDGET_INFO(id.wdropl_datatype, /DROPLIST_SELECT)
	(*id.ptr_st).typedata = id.arr_types[type_select]
	bytes_per_data        = id.arr_bytestype[type_select]

	;-------------------------------------------------------------------------------------------------
	id.size_file_out = bytes_per_data * (1L) * (((*id.ptr_st).xsize)*1L*((*id.ptr_st).ysize)*1L* $
		((*id.ptr_st).nslices)*((*id.ptr_st).nframes)) $
		 + offset
	;-------------------------------------------------------------------------------------------------

	str_tab = STRING(9b)
	WIDGET_CONTROL, id.wlabel_info[1], SET_VALUE= 'Bytes to read through: ' + STRTRIM(id.size_file_out,2) + '   '

	IF id.opt_file THEN BEGIN
		WIDGET_CONTROL, id.wlabel_info[0], SET_VALUE= 'File size (bytes): ' + STRTRIM(id.size_file_bytes,2)+ '   '
    	WIDGET_CONTROL, id.wbttn_ok, SENSITIVE = id.size_file_out LE id.size_file_bytes
    ENDIF

	WIDGET_CONTROL, id.wbase_main, SET_UVALUE= id

	PRINT, 'Event..'
	RETURN, 1


END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_ReadRaw_event, ev

opt_catcherror = 0
;*****************************************
IF LMGR(/RUNTIME) OR LMGR(/VM) OR (opt_catcherror EQ 1) THEN BEGIN
	CATCH, theError
	IF theError NE 0 THEN BEGIN
	  CATCH, /Cancel
	  void = DIALOG_MESSAGE(!ERROR_STATE.MSG, /ERROR)
	  RETURN
	ENDIF
ENDIF
;*****************************************


WIDGET_CONTROL, ev.top, GET_UVALUE=id
WIDGET_CONTROL, ev.id,  GET_UVALUE=uval


IF N_ELEMENTS(uval) GT 0 THEN BEGIN
CASE uval OF
'VALUES': BEGIN
	ok = fwidget_ReadRaw_Event(id)
	WIDGET_CONTROL, id.wbase_main, SET_UVALUE= id
END
'EXIT2':BEGIN
	*id.ptr_st = -1
	Interface_datain_generic_destroy, ev.top
END
'EXIT1':BEGIN
	ok = fwidget_ReadRaw_Event(id)
	WIDGET_CONTROL, id.wbase_main, SET_UVALUE= id
	Interface_datain_generic_destroy, ev.top
END
ELSE:
ENDCASE
ENDIF

END

;**************************************************************************************************
;**************************************************************************************************

PRO Interface_ReadRaw_destroy, ev

	WIDGET_CONTROL, ev, GET_UVALUE=id
	WIDGET_CONTROL,  id.wbase_main,  /DESTROY

END
;**************************************************************************************************
;**************************************************************************************************

FUNCTION Interface_ReadRaw, GROUP_LEADER=group_leader,$
	POSITION=position, FILE=file,$
    ST_DATA=st_data

tlb_frame_attr=1


IF N_ELEMENTS(group_leader) NE 0 THEN BEGIN
    base_floating = WIDGET_BASE(GROUP_LEADER=group_leader, FLOATING=1, UVALUE='INFO_FLOATING', $
       TITLE='Read binary file', TLB_FRAME_ATTR=tlb_frame_attr, MAP=opt_map, COLUMN=1)
ENDIF ELSE BEGIN
    base_floating = WIDGET_BASE(UVALUE='INFO_FLOATING', $
       TITLE='Read binary file', TLB_FRAME_ATTR=tlb_frame_attr, MAP=opt_map, COLUMN=1)
ENDELSE
;-------------------------------
opt_file = 0l
IF N_ELEMENTS(file) NE 0 THEN BEGIN
	IF FILE_TEST(file[0]) EQ 1 THEN BEGIN
		opt_file = 1
		size_file_bytes = (FILE_INFO(file[0])).size
	ENDIF
ENDIF
IF opt_file EQ 0 THEN size_file_bytes = -1
;-------------------------------

base_0 = WIDGET_BASE(base_floating, COLUMN=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)

base_1 = WIDGET_BASE(base_0, ROW=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT, /FRAME)
	base_1a = WIDGET_BASE(base_1, COLUMN=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_1a1 = WIDGET_BASE(base_1a, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT)
		base_1a2 = WIDGET_BASE(base_1a, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT)
		base_1a3 = WIDGET_BASE(base_1a, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT)
		base_1a4 = WIDGET_BASE(base_1a, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT)
	base_1b = WIDGET_BASE(base_1, COLUMN=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_1b1 = WIDGET_BASE(base_1b, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT, /NONEXCLUSIVE)
		base_1b2 = WIDGET_BASE(base_1b, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT, /NONEXCLUSIVE)
		base_1b3 = WIDGET_BASE(base_1b, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT, /NONEXCLUSIVE)
		base_1b4 = WIDGET_BASE(base_1b, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT, /NONEXCLUSIVE)
	base_rell = WIDGET_BASE(base_1, COLUMN=1, XSIZE=10)
	base_1c = WIDGET_BASE(base_1, COLUMN=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_1c1 = WIDGET_BASE(base_1c, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT)
		base_1c2 = WIDGET_BASE(base_1c, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT)
		base_1c3 = WIDGET_BASE(base_1c, ROW=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_1c4 = WIDGET_BASE(base_1c, ROW=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
		base_1c5 = WIDGET_BASE(base_1c, ROW=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)

base_2 = WIDGET_BASE(base_0, ROW=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
	base_2a = WIDGET_BASE(base_2, ROW=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
	base_2b = WIDGET_BASE(base_2, ROW=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)
base_3 = WIDGET_BASE(base_0, ROW=1, /ALIGN_LEFT, /BASE_ALIGN_LEFT)

base_4 = WIDGET_BASE(base_0, /ROW, /ALIGN_CENTER, /BASE_ALIGN_CENTER)


label_xsize   = WIDGET_LABEL(base_1a1, VALUE='Pixels (x)')
text_xsize    = WIDGET_TEXT(base_1a1, /EDITABLE, VALUE=STRTRIM(STRING(st_data.xsize),2),   XSIZE=4, $
	UVALUE='VALUES', UNAME='XSIZE')
label_ysize   = WIDGET_LABEL(base_1a2, VALUE='Pixels (y)')
text_ysize    = WIDGET_TEXT(base_1A2, /EDITABLE, VALUE=STRTRIM(STRING(st_data.ysize),2),   XSIZE=4,$
	UVALUE='VALUES', UNAME='YSIZE')

label_nslices = WIDGET_LABEL(base_1a3, VALUE='Slices (z)')
text_nslices  = WIDGET_TEXT(base_1a3, /EDITABLE, VALUE=STRTRIM(STRING(st_data.nslices),2), XSIZE=4,$
	UVALUE='VALUES', UNAME='SLICES')
label_nframes = WIDGET_LABEL(base_1a4, VALUE='Frames')
text_nframes  = WIDGET_TEXT(base_1a4, /EDITABLE, VALUE=STRTRIM(STRING(st_data.nframes),2), XSIZE=4,$
	UVALUE='VALUES', UNAME='FRAMES')


bttn_invX = WIDGET_BUTTON(base_1b1, VALUE='Invert', UVALUE='VALUES')
bttn_invY = WIDGET_BUTTON(base_1b2, VALUE='Invert', UVALUE='VALUES')
bttn_invZ = WIDGET_BUTTON(base_1b3, VALUE='Invert', UVALUE='VALUES')
bttn_invT = WIDGET_BUTTON(base_1b4, VALUE='Invert', UVALUE='VALUES')

arr_types =     [1,2,3,14, 4, 5,12,13,15L]
arr_bytestype = [1,2,4, 8, 4, 8, 2, 4, 8l] ; bytes per data
strarr_types = [$
	'Byte (8 bits)', $
	'Integer (16 bits)',$
	'Integer (32 bits)',$
	'Integer (64 bits)',$
	'Float (32 bits)', $
	'Float (64 bits)', $
	'Unsigned Int (16 bits)',$
	'Unsigned Int (32 bits)',$
	'Unsigned Int (64 bits)']

dropl_datatype  = WIDGET_DROPLIST(base_1c1, VALUE=strarr_types, TITLE='Data type', $
	UVALUE='VALUES', UNAME='DATA_TYPE')

base_1c2exc = WIDGET_BASE(base_1c2, ROW=1, /ALIGN_RIGHT, /BASE_ALIGN_LEFT, /EXCLUSIVE)
bttn_bigendian = WIDGET_BUTTON(base_1c2exc, VALUE = 'Big endian', UVALUE='VALUES')
bttn_littleendian = WIDGET_BUTTON(base_1c2exc, VALUE = 'Little endian', UVALUE='VALUES')

label_offset = WIDGET_LABEL(base_1c3, VALUE='Offset (bytes)')
text_offset = WIDGET_TEXT(base_1c3, /EDITABLE, VALUE=STRTRIM(STRING(st_data.offset),2), XSIZE=4,$
	UVALUE='VALUES', UNAME='FRAMES')


str_1 = '                            '
label_info1 = WIDGET_LABEL(base_1c4, VALUE='File size (bytes):' + str_1 ,/DYNAMIC_RESIZE)
label_info2 = WIDGET_LABEL(base_1c5, VALUE='Bytes to read through:' + str_1, /DYNAMIC_RESIZE)


label_frame_period = WIDGET_LABEL(base_2a, VALUE='Frame Interval (s)')
text_frame_period  = WIDGET_TEXT(base_2a, /EDITABLE, VALUE=STRTRIM(STRING(st_data.frame_period, FORMAT='(F8.3)'),2),$
	XSIZE=6, UVALUE='VALUES', UNAME='FRAMETIME')

label_injection = WIDGET_LABEL(base_2b, VALUE='Injection frame')
text_injection  = WIDGET_TEXT(base_2b, /EDITABLE, VALUE=STRTRIM(STRING(st_data.frame_injection+1),2),$
	XSIZE=4, UVALUE='VALUES', UNAME='INJECTIONFRAME')

bttn_ok  = WIDGET_BUTTON(base_4, VALUE=' OK     ', UVALUE='EXIT1')
bttn_cancel = WIDGET_BUTTON(base_4, VALUE=' Cancel ', UVALUE='EXIT2')

IF N_ELEMENTS(position) EQ 1 AND N_ELEMENTS(group_leader) NE 0 THEN BEGIN
	pos = GetWidget_RelativePosition(base_floating, ID_PARENT=group_leader, POSITION=position)
	WIDGET_CONTROL, base_floating, XOFFSET=pos[0], YOFFSET=pos[1]
ENDIF ELSE BEGIN
	 pos = GetWidget_PositionInScreen(base_floating,CENTER=1)
ENDELSE

st = { $
	wbase_main    : base_floating,$
	wbttn_ok   		: bttn_ok,$
	wtext_xsize   : text_xsize,$
	wtext_ysize   : text_ysize,$
	wtext_nframes : text_nframes,$
	wtext_nslices : text_nslices,$
	wbttn_inv     : [bttn_invX,bttn_invY,bttn_invZ, bttn_invT],$
	wdropl_datatype:dropl_datatype,$
	wbttn_bigendian : bttn_bigendian,$
	wbttn_littlendian : bttn_littleendian,$
	wtext_offset  : text_offset,$
	wtext_frame_period : text_frame_period,$
	wtext_injection : text_injection,$
	wlabel_info   : [label_info1,label_info2],$

	arr_types       : arr_types,$
	arr_bytestype   : arr_bytestype,$
	opt_file        : opt_file,$
	size_file_bytes : size_file_bytes,$
	size_file_out   : 0l,$
	ptr_st        : PTR_NEW(/ALLOCATE_HEAP),$
	last : 0l }

FOR i=0l, 3 DO BEGIN
	WIDGET_CONTROL, st.wbttn_inv[i], SET_BUTTON   = st_data.invarr[i]
ENDFOR

WIDGET_CONTROL, st.wbttn_littlendian, SET_BUTTON  = st_data.littleendian
WIDGET_CONTROL, st.wbttn_bigendian,   SET_BUTTON  = st_data.littleendian EQ 0

pos_droplist = (WHERE(st_data.Typedata EQ arr_types,ct))[0]
IF ct NE 1 THEN RETURN, -1
WIDGET_CONTROL, st.wdropl_datatype, SET_DROPLIST_SELECT=pos_droplist

*st.ptr_st = st_data

WIDGET_CONTROL, base_floating, SET_UVALUE=st
WIDGET_CONTROL, base_floating, /REALIZE


WIDGET_CONTROL, st.wtext_xsize, SEND_EVENT={ID:0L, TOP:0L, HANDLER:0L}
;ok = fwidget_ReadRaw_Event(id) ; First execution of events..

XMANAGER,  'Interface_ReadRaw',  base_floating, NO_BLOCK=1, MODAL=1, $
    EVENT_HANDLER='Interface_ReadRaw_event', CLEANUP='Interface_ReadRaw_Destroy'

PRINT, 'OK'

result = *st.ptr_st
PTR_FREE, st.ptr_st

RETURN, result

END


;**************************************************************************************************
;**************************************************************************************************


PRO Test_Interface_ReadRaw

st_data = {xsize:128l, ysize:128l, nslices:4l, nframes:40l, offset:0, $
frame_period:40.0, frame_injection:11, typedata:2, $
invarr : [0,1,0,0l], littleendian:0}


str_file = 'F:\Experiments_DCE-MRI\dMRI_phantoms\' + 'Phantom_Hoffmann_v1.raw'

std = Interface_ReadRaw(GROUP_LEADER=group_leader,FILE=str_file, POSITION=[20,40], ST_DATA=st_data)

PRINT, std
PRINT, ''
PRINT, st_data
PRINT, ''

IF SIZE(std, /TNAME) NE 'STRUCT' THEN RETURN

dimensions = [std.xsize, std.ysize, std.nslices, std.nframes]
data = Read_raw_v2(FILE=str_file, OFFSET=std.offset, DIMENSIONS=dimensions, SWAP=std.littleendian, $
		TYPE=std.typedata, REVERSE_ARR=std.invarr, ERROR_ID = error_id, ERROR_STR = error_str,$
		MATCH_SIZE=0)

HELP, data

VIEWG, data[*,*,1,16]


END

;<-->
;<+>

FUNCTION  Locate_FirstNoValid_Widget, array_widgets

n_pos = N_ELEMENTS(array_widgets)

FOR i=0, n_pos-1 DO BEGIN
	valid_id = WIDGET_INFO(array_widgets[i], /VALID_ID)
	IF valid_id EQ 0 THEN RETURN, i
ENDFOR

RETURN, -1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<++>

PRO main

	only_one_instance = 0

	n_instance =  XRegistered('Interface_DCEMRI')
	IF only_one_instance AND n_instance GT 0 THEN RETURN

	filename_ico  = 'logohead_48x48.ico'
	file_bitmap_1 = '.\Icons\' + filename_ico
	file_bitmap_2 = '.\' + filename_ico
	file_bitmap_3 = '..\Icons\' + filename_ico
	IF FILE_TEST(file_bitmap_1) THEN image_bitmap = file_bitmap_1
	IF FILE_TEST(file_bitmap_2) THEN image_bitmap = file_bitmap_2
	IF FILE_TEST(file_bitmap_3) THEN image_bitmap = file_bitmap_3

	;im = read_png('.\data\logo_BIT_head_32x32.png')
	;im = change_true(im)

	DCEMRI_WelcomeWidget, OK=ok
	IF ok EQ 0 THEN RETURN

	Interface_DCEMRI, IMAGE_BITMAP=image_bitmap

END

;<-->
;<+>

FUNCTION Model_Cp_Biexponential, DOSE=dose, ARR_TIME=arr_time, PARAMS=params, STRUCT=struct, HAEMATOCRIT=haematocrit


	IF N_ELEMENTS(dose) EQ 0 THEN dose_temp = 1d ELSE dose_temp = dose[0]
	;a1 = (kg/litro)
	;a2 = (kg/litro)
	;m1 = min-1 o seg-1
	;m2 = min-1 o seg-1
	; arr_time (en segundos o minutos de acuerdo a m1,m2)
	; dose = mmol/kg animal

	IF N_ELEMENTS(struct) EQ 0 THEN BEGIN
		a1 = params[0]
		a2 = params[1]
		m1 = params[2]
		m2 = params[3]
	ENDIF ELSE BEGIN
		a1 = struct.a1[0]
		a2 = struct.a2[0]
		m1 = struct.m1[0]
		m2 = struct.m2[0]
	ENDELSE

	arr_C_p = Dose_temp*(a1*exp(-m1*arr_time) + a2*exp(-m2*arr_time)); Ecuacion 1 (Tofts ISMRM) tiempo en minutos

	;arr_C_p[0] = 0

	IF N_ELEMENTS(haematocrit) NE 0 THEN BEGIN
		arr_C_p = arr_C_p/((1.0d)-haematocrit[0])
	ENDIF

	RETURN, arr_C_p
	; resultado en t(0) = dose*(a1+a2) mmol/kg*(kg/litro) = mmol/litro o mM

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Model_Cp_Orton2, ARR_TIME_MIN=arr_time_min, PARAMS=params, STRUCT=struct, HAEMATOCRIT=haematocrit

; Modelos 2 para la AIF de:
; Orton et al., Computationally efficient vascular input function models for quantitative
; kinetic modelling using DCE-MRI. Phys Med Biol 2008, 53(5):1225-1239.

; (Función triexponencial)

; numeros del articulo:
;ab = 344    mM/min;
;ag = 1.24   min-1;
;ub = 20.2   min-1,
;ug = 0.172  min-1

IF N_ELEMENTS(struct) EQ 0 THEN BEGIN
	ab = params[0]
	ag = params[1]
	ub = params[2]
	ug = params[3]
ENDIF ELSE BEGIN
	ab = struct.ab[0]
	ag = struct.ag[0]
	ub = struct.ub[0]
	ug = struct.ug[0]
ENDELSE


A_b = ab - ab*ag/(ub-ug)

A_g = ab*ag/((ub-ug)^2)

arr_exp_ub = EXP(-ub*arr_time_min)

Arr_Cp = A_b*arr_time_min*arr_exp_ub + A_g*(EXP(-ug*arr_time_min) - arr_exp_ub)

RETURN, arr_Cp


END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Model_Cp_Parker, ARR_TIME_MIN=arr_time_min, HAEMATOCRIT=haematocrit, PARAMS=params, STRUCT=struct

; modelo de Parker para la concentracion en plasma

; OPCION 1:
;------------------------------------------------------------------------------------
; G. J. M. Parker et al.
; "Experimentally-derived functional form for a population-averaged high-temporal-resolution
;  arterial input function for dynamic contrast-enhanced MRI"
;  Magnetic Resonance in Medicine, vol. 56, (5), pp. 993-1000, 2006
;------------------------------------------------------------------------------------
; Tambien utilizado en Guo et al
;(Guo et al, 2009)
;-----------------------------------------------------------------------------------------
; J. Y. Guo, et al. Dynamic contrast-enhanced magnetic resonance imaging parameters independent
; of baseline T-10 values”, Magnetic Resonance Imaging, vol. 27, no. 9, pp. 1208-1215, Nov, 2009.
;-----------------------------------------------------------------------------------------


; La expresión de Guo et al, 2009 tiene distinto significado para A1, A2

; HTC= hematocrito

;alpha = 1.05d     ;Unidades: mmol
;beta  = 0.1685d   ;Unidades: Min-1
;s   = 38.078d     ;Unidades: Min-1
;tau = 0.483d      ;Unidades: Min
;tau1 = 0.17046d   ;Unidades: Min
;tau2 = 0.365d     ;Unidades: Min
;sigma1 = 0.0563d  ;Unidades: Min
;sigma2 = 0.132d   ;Unidades: Min
;A1 = 0.809d       ;Unidades: mmol*min
;A2 = 0.330d       ;Unidades: mmol*min

arr_t_min = DOUBLE(arr_time_min)

IF N_ELEMENTS(struct) EQ 0 THEN BEGIN
	IF N_ELEMENTS(params) EQ 0 THEN BEGIN
	 ;(Parker et al. 2006)

		alpha = 1.05d     ; Unidades: mmol
		beta  = 0.1685d   ; Unidades: Min-1

		s   = 38.078d     ;Unidades: Min-1
		tau = 0.483d      ;Unidades: Min

		tau1 = 0.17046d   ;Unidades: Min
		tau2 = 0.365d     ;Unidades: Min

		sigma1 = 0.0563d  ;Unidades: Min
		sigma2 = 0.132d   ;Unidades: Min

		A1 = 0.809d       ;Unidades: mmol*min
		A2 = 0.330d       ;Unidades: mmol*min
	ENDIF ELSE BEGIN
		alpha  = params[0]    ; Unidades: mmol
		beta   = params[1]    ; Unidades: Min-1
		s      = params[2]    ; Unidades: Min-1
		tau    = params[3]    ; Unidades: Min
		tau1   = params[4]    ; Unidades: Min
		tau2   = params[5]    ; Unidades: Min
		sigma1 = params[6]    ; Unidades: Min
		sigma2 = params[7]    ; Unidades: Min
		A1     = params[8]    ; Unidades: mmol*min
		A2     = params[9]    ; Unidades: mmol*min
	ENDELSE
ENDIF ELSE BEGIN
		alpha  = struct.alpha[0]  ; Unidades: mmol
		beta   = struct.beta[0]   ; Unidades: Min-1
		s      = struct.s[0]      ; Unidades: Min-1
		tau    = struct.tau[0]    ; Unidades: Min
		tau1   = struct.t1[0]     ; Unidades: Min
		tau2   = struct.t2[0]     ; Unidades: Min
		sigma1 = struct.sigma1[0] ; Unidades: Min
		sigma2 = struct.sigma2[0] ; Unidades: Min
		A1     = struct.a1[0]     ; Unidades: mmol*min
		A2     = struct.a2[0]     ; Unidades: mmol*min
ENDELSE


arr_C_p = alpha*EXP(-beta*arr_t_min)/(1+ EXP(-s*(arr_t_min-tau))) + $   ; caída lenta (equiv al biexponencial de tofts)
		 (A1/SQRT(2*!DPI)/sigma1)*EXP(-(arr_t_min-tau1)^2/(2.0d*(sigma1^2))) + $ ; segundo lobulo (recirculación)
		 (A2/SQRT(2*!DPI)/sigma2)*EXP(-(arr_t_min-tau2)^2/(2.0d*(sigma2^2)))     ; primer lóbulo (subida principal)


IF N_ELEMENTS(haematocrit) NE 0 THEN BEGIN
	arr_C_p = arr_C_p/((1.0d)-haematocrit[0])
ENDIF

RETURN, arr_C_p


END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Model_Cp_Schabel, ARR_TIME_MIN=arr_time_min, OPTION=option, DELTA_T_MIN=delta_t_min, $
			ST_PARAM=st_param, $
			ST_IN=st_in, SEED=seed, $
			PARAMS_10_OUT=params_10_out,$
			PARAMS_10_IN =params_10_in, $
			A0=A0


; La entrada puede ser St_in, PARAMS_10_IN  o bien OPTION
;

; ARR_TIME_SEG= arr_time_Seg
; DELTA_T_SEG = comienzo de la subida, en segundos


; Función AIF (Arterial input function) de forma gamma-variate + sigmoid
; Según (Schabel et al. 2008, Apéndice A)
;---------------------------------------------------------------------------------
;M. C. Schabel, and D. L. Parker,
;Uncertainty and bias in contrast concentration measurements using spoiled gradient echo pulse sequences,
;Physics in Medicine and Biology, vol. 53, no. 9, pp. 2345-2373, May 7, 2008.
;---------------------------------------------------------------------------------

; Ver  también
;(Fluckiger et. al., 2009) ecuación 3
;---------------------------------------------------------------------------------
; J. U. Fluckiger, M. C. Schabel, and E. V. R. DiBella,
; Model-Based Blind Estimation of Kinetic Parameters in Dynamic Contrast Enhanced (DCE)-MRI
; Magnetic Resonance in Medicine, vol. 62, no. 6, pp. 1477-1486, Dec, 2009.
;---------------------------------------------------------------------------------

; Y también descrito en (Schabel et al, 2010) (parte 1)
;---------------------------------------------------------------------------------
; M. C. Schabel, J. U. Fluckiger, and E. V. R. DiBella,
; "A model-constrained Monte Carlo method for blind arterial input function estimation
;  in dynamic contrast-enhanced MRI Simulations,"
; Physics in Medicine and Biology, vol. 55, (16), pp. 4783-4806, 2010.
;---------------------------------------------------------------------------------

n_points = N_ELEMENTS(arr_time_min)

arr_t_min = arr_time_min
IF N_ELEMENTS(delta_t_min) EQ 0 THEN BEGIN
	delta_t_min   = 0d
ENDIF ELSE BEGIN
	delta_t_min  = delta_t_min[0]
ENDELSE

; cuidado con el orden, la sigmoide va al final en este codigo, pero en (Schabel et al. 2008)
; y (Fluckinger et al, 2009) iría al principio

opt_10_params = 0l

IF N_ELEMENTS(st_in) GT 0 THEN BEGIN
	; se proporcionan los datos en st_in
	A       = st_in.A
	alpha   = st_in.alpha
	tau     = st_in.tau
	delay_t = st_in.delay_t
	delta_t = st_in.delta_t + delta_t_min
ENDIF ELSE IF N_ELEMENTS(params_10_in) EQ 10 THEN BEGIN ; se dan solamente 10 parámetros propuestos por Schabel et al. (2010)

	;params  = [ A[1]/A[0], A[2]/A[0], A[3]/A[0],  alpha[0],$
	;		    tau[0], tau[1], delay_t[0], delta_t[0], delta_t[1]-delta_t[0], delta_t[2]-delta_t[0]]

	A       = [1.0, params_10_in[0:2]]
	alpha   = REPLICATE(params_10_in[3],4)
	tau     = [params_10_in[4], REPLICATE(params_10_in[5],3)]
	delay_t =  params_10_in[6]
	delta_t = [params_10_in[7], params_10_in[8]+params_10_in[7], $
		params_10_in[9]+params_10_in[7], params_10_in[9]+params_10_in[7]]

ENDIF ELSE BEGIN

;  hay que utilizar una opción por defecto para generar st_param

	IF N_ELEMENTS(option) EQ 0 THEN RETURN, -1 ELSE opt=FIX(option[0])

	CASE opt OF

	; (Schabel et al. 2008)
	1 : BEGIN ; (Schabel et al. 2008) three normalized gamma-variates and 1 sigmoid
		A     = [5.8589,  0.9444,  0.4888, 0.8152d]
		alpha = [2.5393,  7.9461,  7.9461, 7.9461d]
		tau   = [0.04286, 0.06873, 0.1400, 0.14d] ; en minutos
		delay_t = 9.6319d ; en minutos
		delta_t = REPLICATE(delta_t_min[0],4)*1d
	END
	2 : BEGIN ;(Schabel et al. 2010) Cpop (population averages)
		A       = [6.0, 1.1208, 0.3024, 0.7164d]  ; mM (milimolar)
		alpha   = [2.92, 2.92, 2.92, 2.92d]       ; adimensional (alpha_1=alpha_2=alpha_3=alpha_4
		tau     = [0.0442, 0.1430, 0.1430, 0.1430d]; min (tau_2=tau_3=tau_4)
		delay_t = 7.8940      ; en minutos, T mayuscula en el articulo
		delta_t = [0, 0.2227, 0.6083, 0.6083]+delta_t_min  ; en minutos,  delta_4=delta_3
		; los diez parametros que utiliza "Model_Schabel_Cp_v2", para chequeo
		opt_10_params = 1l

	END
	3 : BEGIN ;(Schabel et al. 2010) Cmod (modelada)
		A       = [6.0, 1.6200, 0.7164,1.2732]  ; mM (milimolar)
		alpha   = [2.92, 2.92, 2.92, 2.92d]       ; adimensional (alpha_1=alpha_2=alpha_3=alpha_4
		tau     = [0.0633, 0.0954, 0.0954d, 0.0954d]; min (tau_2=tau_3=tau_4)
		delay_t = 14.7120d      ; en minutos, T mayuscula en el articulo
		delta_t = [0.0, 0.300, 0.6083, 0.6083d]+delta_t_min  ; en minutos,  delta_4=delta_3
		opt_10_params = 1l
	END
	4 : BEGIN ; (Fluckiger et. al., 2009)		; two normalized gamma-variates and 1 sigmoid
		A       = [7.8, 1.2, 1.9d]             ;
		alpha   = [1.26, 1.65, 44.0d]
		tau     = [0.00663, 0.0055, 0.018d] ; en minutos
		delay_t = 15.0d    ;en minutos
		delta_t = REPLICATE(delta_t_min[0],3)*1d  ; 0.7
	END
	5 : BEGIN ; (Fluckiger et. al., 2009)		; variación
		A       = [4, 1.5, 2.0d]
		alpha   = [1.1, 1.5, 40.0d]
		tau     = [0.009, 0.01, 0.03] ; en minutos
		delay_t = 20.0d
		delta_t = REPLICATE(delta_t_min[0],3)*1d  ; 0.5
	END
	6 : BEGIN
		IF N_ELEMENTS(seed) EQ 0 THEN seed=0
		; de prueba...(variacion aleatoria respecto a 2. (pagina 4789, schabel et al, 2010)
		A = DBLARR(4)
		A[0] = 1.0          ; Se deja fijo el valor máximo del lóbulo principal
		A[1] = A[0]*RANDOMU(seed,1)*0.2
		A[2] = A[1]*RANDOMU(seed,1)*0.5
		A[3] = A[0]*RANDOMU(seed,1)*0.5

		alpha = REPLICATE(((RANDOMU(seed,1))[0]+1)*4, 4) ; entre 1 y 5
		delta_t = DBLARR(4)
		delta_t[1] = (RANDOMU(seed,1))[0] ; entre 0 y 1
		delta_t[2] = delta_t[1] + (RANDOMU(seed,1))[0] ; entre delta_t[2] + [0,1]
		delta_t[3] = delta_t[2] + (RANDOMU(seed,1))[0] ; entre delta_t[2] + [0,1]
		delta_t+=delta_t_min
		tau = DBLARR(4)
		tau[0] = RANDOMU(seed,1)*0.2
		tau[[1,2,3]] = REPLICATE(RANDOMU(seed,1)*1.2, 3)
		delay_t = RANDOMU(seed,1)*10 + tau[0]

		;PRINT, 'A', A
		;PRINT, 'alpha', alpha
		;PRINT, 'tau', tau
		;PRINT, 'delay_t', deLay_t
		;PRINT, 'delta_t', delta_t
		;PRINT, ''
		opt_10_params = 1
	END
	ENDCASE
ENDELSE


n_sigmas = N_ELEMENTS(A)-1
pos_sig  = n_sigmas

data1 = A[pos_sig]*Sigmoid_curve(arr_t_min-delta_t[pos_sig], DELAY_T=delay_t, ALPHA=alpha[pos_sig], TAU=tau[pos_sig]) ; el sigmoide

data_temp = DBLARR(n_points)

FOR i=0l, n_sigmas-1 DO BEGIN
	n = i
	data_temp+=A[n]*NormalizedGamma_variate(arr_t_min-delta_t[n], ALPHA=alpha[n], TAU=tau[n])
ENDFOR

data_aif = data1+data_temp

st_param = {A:A,alpha:alpha,tau:tau, delay_t:delay_t, delta_t:delta_t}
;----------------------------------------------------------------------------------------
IF opt_10_params EQ 1 THEN BEGIN
	; los diez parametros que utiliza "Model_Schabel_Cp_v2", para chequeo
	params_10_out  = [ A[1]/A[0], A[2]/A[0], A[3]/A[0],  alpha[0],$
		tau[0], tau[1], delay_t[0], delta_t[0], delta_t[1]-delta_t[0], delta_t[2]-delta_t[0]]
ENDIF
;----------------------------------------------------------------------------------------

A0 = A[0]

RETURN, data_aif

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Model_Cp_SchabelSimplified, ARR_TIME_MIN=arr_time_min, PARAMS=params, STRUCT=struct, HAEMATOCRIT=haematocrit

; En realidad ajusta una variación de Schabel con una gamma y una sigmoide (schabel simplified)

IF N_ELEMENTS(struct) EQ 0 THEN BEGIN
	A1      = params[0]
	A2      = params[1]
	alpha   = params[2]
	tau1    = params[3]
	tau2    = params[4]
	delay_t = params[5]
ENDIF ELSE BEGIN
	A1      = struct.a1[0]
	A2      = struct.a2[0]
	alpha   = struct.alpha[0]
	tau1    = struct.tau1[0]
	tau2    = struct.tau2[0]
	delay_t = struct.delay[0]
ENDELSE

data_sigmoid   = A2*Sigmoid_curve(arr_time_min, DELAY_T=delay_t, ALPHA=alpha, TAU=tau2) ; el sigmoide
data_normgamma = A1*NormalizedGamma_variate(arr_time_min, ALPHA=alpha, TAU=tau1)

RETURN, data_sigmoid + data_normgamma

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Model_Ct_Tofts, DOSE=dose, ARR_TIME=arr_time, KTRANS=ktrans, KEP=Kep,  VP=vp, $
	TAU_SHIFT=tau_shift, ARR_AUX=arr_aux, PARAMS=params, MODEL_AIF=model_aif, IDS_INPUT=ids_input

; Modelo de Tofts apra concentracion en tejido (Ct) con concentración en plasma (Cp) biexponencial

;TAU_SHIFT  -> desfase adicional de la biexponencial Cp (Orton et al, 2007)
;VP -> Si es distinto de cero, realiza el modelo de Tofts modificado con término Vp*Cp


IF N_ELEMENTS(ids_input) GT 0 THEN BEGIN
	param_ids     = ['KTRANS', 'KEP', 'VP', 'TAU_SHIFT']
	npos_needed   = WHERE([1,1,0,0l] EQ 1)
	params_exist  = INTARR(4)
	params_f      = DBLARR(4)
	n_input       = N_ELEMENTS(params)
	FOR i=0l, n_input-1 DO BEGIN
		pos = (WHERE(ids_input[i] EQ param_ids, cpos))[0]
		IF cpos EQ 1 THEN BEGIN
			params_f[pos] = params[i]
			params_exist[pos]+=1
		ENDIF
	ENDFOR
	IF MIN(params_exist[npos_needed]) NE 1 THEN RETURN, -1
	ktrans    = params_f[0]
	kep       = params_f[1]
	vp        = params_f[2]
	tau_shift = params_f[3]
ENDIF ELSE BEGIN
	;------------------------------------------------
	IF N_ELEMENTS(params) GT 0 THEN BEGIN ; to integrate in an interface with another modeling functions
		ktrans    = params[0]
		kep       = params[1]
		vp        = params[2]
		tau_shift = params[3]
	ENDIF
ENDELSE
IF N_ELEMENTS(vp) EQ 0 THEN volume_plasma= 0d ELSE volume_plasma = vp[0]
IF N_ELEMENTS(tau_shift) EQ 0 THEN tau_t=0d ELSE tau_t = tau_shift[0]

; modelo de Tofts
;-----------------------------------------------
a1        = model_aif[0] ; (kg/litro)
a2        = model_aif[1] ; (kg/litro)
m1        = model_aif[2] ; s-1 o min-1
m2        = model_aif[3] ; s-1 o min-1
;-----------------------------------------------

; modelo de Tofts


;arr_C_t = D*KTRANS*(a1*(exp(-kep*arr_t_min)-exp(-m1*arr_t_min))/(m1-kep) + $
;					a2*(exp(-kep*arr_t_min)-exp(-m2*arr_t_min))/(m2-kep))

arr_t_shift = arr_time-tau_t

arr_Ct =  Dose*KTRANS*(a1*(exp(-kep*arr_t_shift)-exp(-m1*arr_t_shift))/(m1-kep) + $
					     a2*(exp(-kep*arr_t_shift)-exp(-m2*arr_t_shift))/(m2-kep))

IF volume_plasma GT 0 THEN BEGIN ;Modified Tofts
	arr_Cp = Model_Cp_Biexponential(DOSE=dose, ARR_TIME=arr_t_shift, ARR_AUX=arr_aux)
	arr_Cp[0]=0
	arr_Ct = arr_Ct + volume_plasma[0]*arr_Cp
ENDIF

RETURN, arr_Ct



END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Model_Ct_Tofts_convol, ARR_TIME=arr_time, ARR_CP=arr_cp, $
	KTRANS=ktrans, KEP=kep, VP=vp, NORMALIZE=normalize, METHOD=method, PARAMS=params, IDS_INPUT=ids_input

; Tofts model as a convolution, provided general array of Cp (arr_cp)
; El resultado se aproxima más al analítico cuando el intervalo temporal de frames es pequeño
;
; Method:
; 1 - Convolution
; 2 - Convolution matrix
; 3 - Truncated convolution matrix

IF N_ELEMENTS(tau_shift) EQ 0 THEN tau_t = 0d ELSE tau_t = tau_shift[0]

n_times = N_ELEMENTS(arr_time)
IF N_ELEMENTS(arr_cp) NE n_times THEN RETURN, -1


IF N_ELEMENTS(ids_input) GT 0 THEN BEGIN
	param_ids = ['KTRANS', 'KEP', 'VP']
	npos_needed   = WHERE([1,1,0] EQ 1)
	params_exist  = INTARR(3)
	params_f      = DBLARR(3)
	n_input       = N_ELEMENTS(params)
	FOR i=0l, n_input-1 DO BEGIN
		pos = (WHERE(ids_input[i] EQ param_ids, cpos))[0]
		IF cpos EQ 1 THEN BEGIN
			params_f[pos] = params[i]
			params_exist[pos]+=1
		ENDIF
	ENDFOR
	IF MIN(params_exist[npos_needed]) NE 1 THEN RETURN, -1
	ktrans   = params_f[0]
	kep      = params_f[1]
	vp       = params_f[2]
ENDIF ELSE BEGIN
	;------------------------------------------------
	IF N_ELEMENTS(params) GT 0 THEN BEGIN ; to integrate in an interface with another modeling functions
		ktrans  = params[0]
		kep     = params[1]
		vp      = params[2]
	ENDIF
ENDELSE
IF N_ELEMENTS(vp) EQ 0 THEN vp=0d


kernel    = ktrans[0]*EXP(-kep[0]*arr_time)
norm_time = (ABS(arr_time[0]-arr_time[N_ELEMENTS(arr_time)-1]))/(N_ELEMENTS(arr_time)-1)

IF N_ELEMENTS(method) EQ 0 THEN opt_method=1 ELSE opt_method=method[0]

CASE opt_method OF
1 : BEGIN ; convolucion

	arr_zeros = FLTARR(n_times)
	arr_Ctt    = CONVOL([arr_zeros, arr_Cp], REVERSE(kernel), NORMALIZE=normalize)

	p1 = n_times - n_times/2
	p2 = n_times + n_times/2 - 1

	IF ((n_times MOD 2) EQ 0) THEN BEGIN
		p1++ & p2++ ; cambio para compatibilidad MATLAB (sumar 1) caso par
	ENDIF; caso par

	IF ((p2-p1+1) LT n_times) THEN p2++ ; caso de arrays impares

	;pos_nz = WHERE(arr_Ct GT 0, n_pos)

	arr_Ct = (arr_Ctt[p1: p2])*norm_time

	IF N_ELEMENTS(vp) GT 0 THEN BEGIN
		arr_Ct = arr_Ct + vp[0]*arr_Cp
	ENDIF
END
2 : BEGIN
	; matriz de convolución

	Hmatrix = ConvolutionMatrix(kernel*norm_time, n_times) ; matriz de convolucion

	DiagM = DIAG_MATRIX(REPLICATE(vp, n_times));
	DiagM = [[DiagM], [DBLARR(n_times, n_times-1)]]

	arr_Ct = (Hmatrix+DiagM)##TRANSPOSE(arr_Cp)
	arr_Ct = arr_Ct[0:n_times-1]
END
3 : BEGIN
	; matriz de convolución pero cortada...
	Hmatrix = ConvolutionMatrix(kernel*norm_time, n_times) ; matriz de convolucion
	Hmatrix = Hmatrix[0:n_times-1,0:n_times-1]
	DiagM = DIAG_MATRIX(REPLICATE(vp, n_times));
	arr_Ct = (Hmatrix+DiagM)##TRANSPOSE(arr_Cp)
	arr_Ct = REFORM(arr_Ct)
END
ELSE : RETURN, -1
ENDCASE

RETURN, arr_Ct

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Model_Ct_Yankeelov_convol, ARR_TIME=arr_time, ARR_CT_RR=arr_ct_rr, $
	KTRANS_T=ktrans_t, KEP_t=kep_t, KTRANS_RR=ktrans_rr, KEP_RR=kep_rr, NORMALIZE=normalize, METHOD=method,$
	PARAMS=params, IDS_INPUT=ids_input

; Tofts model as a convolution, provided general array of Cp (arr_cp)
; El resultado se aproxima más al analítico cuando el intervalo temporal de frames es pequeño
;
; Method:
; 1 - Convolution
; 2 - Convolution matrix
; 3 - Truncated convolution matrix

IF N_ELEMENTS(tau_shift) EQ 0 THEN tau_t = 0d ELSE tau_t = tau_shift[0]

n_times = N_ELEMENTS(arr_time)
IF N_ELEMENTS(arr_ct_rr) NE n_times THEN RETURN, -1


IF N_ELEMENTS(ids_input) GT 0 THEN BEGIN
	param_ids = ['KTRANS', 'KEP']
	npos_needed   = WHERE([1,1l] EQ 1)
	params_exist  = INTARR(2)
	params_f      = DBLARR(2)
	n_input       = N_ELEMENTS(params)
	FOR i=0l, n_input-1 DO BEGIN
		pos = (WHERE(ids_input[i] EQ param_ids, cpos))[0]
		IF cpos EQ 1 THEN BEGIN
			params_f[pos] = params[i]
			params_exist[pos]+=1
		ENDIF
	ENDFOR
	IF MIN(params_exist[npos_needed]) NE 1 THEN RETURN, -1
   	ktrans_t  = params_f[0]
	kep_t     = params_f[1]
ENDIF ELSE BEGIN
	;------------------------------------------------
	IF N_ELEMENTS(params) GT 0 THEN BEGIN ; to integrate in an interface with another modeling functions
		ktrans_t  = params[0]
		kep_t     = params[1]
	ENDIF
	;-----------------------------------------------
ENDELSE

ve_t  = ktrans_t[0]/kep_t[0]
ve_rr = ktrans_rr[0]/kep_rr[0]
vp    = 0d


kernel    = (ktrans_t[0]/ktrans_RR[0])*(ktrans_RR[0]/ve_RR[0] - ktrans_t[0]/ve_t[0])*EXP(-(ktrans_t[0]/ve_t[0])*arr_time)

norm_time = (ABS(arr_time[0]-arr_time[N_ELEMENTS(arr_time)-1]))/(N_ELEMENTS(arr_time)-1)

IF N_ELEMENTS(method) EQ 0 THEN opt_method=1 ELSE opt_method=method[0]

CASE opt_method OF
1 : BEGIN ; convolucion

	arr_zeros = FLTARR(n_times)
	arr_Ct    = CONVOL([arr_zeros, arr_ct_rr], REVERSE(kernel), NORMALIZE=normalize)
	;--------------------------------------------------------------------------
	p1 = n_times - n_times/2
	p2 = n_times + n_times/2 - 1

	IF ((n_times MOD 2) EQ 0) THEN BEGIN
		p1++ & p2++ ; cambio para compatibilidad MATLAB (sumar 1) caso par
	ENDIF; caso par

	IF ((p2-p1+1) LT n_times) THEN p2++ ; caso de arrays impares
	;--------------------------------------------------------------------------
	arr_Ct = (arr_Ct[p1: p2])*norm_time + (ktrans_t[0]/ktrans_RR[0])*arr_ct_rr

END
2 : BEGIN
	; matriz de convolución

	Hmatrix = ConvolutionMatrix(kernel*norm_time, n_times) ; matriz de convolucion

	DiagM = DIAG_MATRIX(REPLICATE(vp, n_times));
	DiagM = [[DiagM], [DBLARR(n_times, n_times-1)]]
	arr_Ct = (Hmatrix+DiagM)##TRANSPOSE(arr_ct_rr)
	arr_Ct = arr_Ct[0:n_times-1] + (ktrans_t[0]/ktrans_RR[0])*arr_ct_rr

END
3 : BEGIN
	; matriz de convolución pero cortada...
	Hmatrix = ConvolutionMatrix(kernel*norm_time, n_times) ; matriz de convolucion
	Hmatrix = Hmatrix[0:n_times-1,0:n_times-1]
	DiagM = DIAG_MATRIX(REPLICATE(vp, n_times));
	arr_Ct = (Hmatrix+DiagM)##TRANSPOSE(arr_ct_rr)
	arr_Ct = REFORM(arr_Ct) + (ktrans_t[0]/ktrans_RR[0])*arr_ct_rr

END
ELSE : RETURN, -1
ENDCASE

RETURN, arr_Ct

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Model_RCE_Hoffmann, ARR_TIME=arr_time, Kel=Kel, Kep=Kep, AVALUE=avalue, TAU_SHIFT=tau_shift,$
	TAU_INFUSION=tau_infusion, BVALUE=bvalue, M1VALUE=m1value, PARAMS=params, IDS_INPUT=ids_input


; Modelo de Hoffmann para contraste relativo (RCE)

;- TAU_SHIFT    - desfase adicional
;- TAU_INFUSION - Modelo con tiempo de infusión

; (formula en Hoffmann et al. 1995, Huang et al. 2004)
; Notación: (Tofts, 1997)


IF N_ELEMENTS(ids_input) GT 0 THEN BEGIN
	param_ids     = ['KEL', 'KEP', 'AH', 'TAU_SHIFT', 'TAU_INFUSION', 'BH', 'M1']
	npos_needed   = WHERE([1,1,1,0l,0,0,0] EQ 1)
	params_exist  = INTARR(7)
	params_f      = DBLARR(7)
	n_input       = N_ELEMENTS(params)
	FOR i=0l, n_input-1 DO BEGIN
		pos = (WHERE(ids_input[i] EQ param_ids, cpos))[0]
		IF cpos EQ 1 THEN BEGIN
			params_f[pos] = params[i]
			params_exist[pos]+=1
		ENDIF
	ENDFOR
	IF MIN(params_exist[npos_needed]) NE 1 THEN RETURN, -1
	kel          = params_f[0]
	kep          = params_f[1]
	Avalue       = params_f[2]
	tau_shift    = params_f[3]
	tau_infusion = params_f[4]
	bvalue       = params_f[5]
	m1value      = params_f[6]
ENDIF ELSE BEGIN
	;------------------------------------------------
	IF N_ELEMENTS(params) GT 0 THEN BEGIN ; to integrate in an interface with another modeling functions
		kel          = params[0]
		kep          = params[1]
		Avalue       = params[2]
		tau_shift    = params[3]
		tau_infusion = params[4]
		bvalue       = params[5]
		m1value      = params[6]
	ENDIF
ENDELSE
;-----------------------------------------------

IF N_ELEMENTS(tau_shift) EQ 0 THEN tau_t = 0d ELSE tau_t = tau_shift[0]
IF N_ELEMENTS(tau_infusion) EQ 0 THEN tau_inf = 0.0 ELSE tau_inf = tau_infusion[0]

IF tau_inf EQ 0 THEN BEGIN

	;expresión en tofts_1997
	arr_rce = Avalue*kep*(exp(-kep*(arr_time-tau_t))-exp(-kel*(arr_time-tau_t)))/(kel-kep) + 1

ENDIF ELSE BEGIN

	rel1 = kep/(kel*(kep-kel))
	rel2 = 1/(kep-kel)
	; expresión en formula de Hoffmann 1995.
	; se ajusta con MPFIT_RCE_HOFFMANN_V3

	arr_tp = (arr_time < tau_inf[0]) ; during infusion (t < tau), tp=t, after infusion, tp=tau

	arr_rce = (Avalue/tau_inf)*(rel1*(exp(kel*arr_tp)-1)*exp(-kel*arr_time)  - $
		         rel2*(exp(kep*arr_tp)-1)*exp(-kep*arr_time) ) + 1
ENDELSE

IF N_ELEMENTS(bvalue) NE 0 THEN BEGIN

	arr_bval = Bvalue*EXP(-m1value[0]*arr_time)
	arr_rce = arr_rce+arr_bval

ENDIF

RETURN, arr_rce

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Model_RCE_Larsson, ARR_TIME=arr_time, PARAMS=params, MODEL_AIF=model_aif, IDS_INPUT=id_input

; (Tofts et al, 1997)


IF N_ELEMENTS(ids_input) GT 0 THEN BEGIN
	param_ids = ['SL', 'KEP']
	npos_needed   = WHERE([1,1l] EQ 1)
	params_exist  = INTARR(2)
	params_f      = DBLARR(2)
	n_input       = N_ELEMENTS(params)
	FOR i=0l, n_input-1 DO BEGIN
		pos = (WHERE(ids_input[i] EQ param_ids, cpos))[0]
		IF cpos EQ 1 THEN BEGIN
			params_f[pos] = params[i]
			params_exist[pos]+=1
		ENDIF
	ENDFOR
	IF MIN(params_exist[npos_needed]) NE 1 THEN RETURN, -1
	kep = params_f[0]
	S   = params_f[1]
ENDIF ELSE BEGIN
	kep = params[0]
	S   = params[1]
ENDELSE

a1  = model_aif[0]
a2  = model_aif[1]
m1  = model_aif[2]
m2  = model_aif[3]

arr_sal = 1.0 + (S/(a1+a2)) * ( a1*(EXP(-kep*arr_time) - EXP(-m1*arr_time))/(m1-kep) + $
					             a2*(EXP(-kep*arr_time) - EXP(-m2*arr_time))/(m2-kep) )

RETURN, arr_sal


END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION MPFIT_CT_TOFTS_BIEXP, X, P

; Ajuste del modelo de Tofts para la función 'MPFITFUN'
; Concentración en tejido (Ct)
;
; Versión con la opción añadida del shift temporal (tau) (Orton et al, 2007)

	n_x = N_ELEMENTS(X)

	;a1 = (kg/litro)
	;a2 = (kg/litro)
	;m1 = min-1
	;m2 = min-1

	D     = X[0]
	a1    = X[1]
	a2    = X[2]
	m1    = X[3]
	m2    = X[4]
	arr_t = X[5:n_x-1]

	KTRANS= P[0]
	ve   =  P[1]
	vp   =  P[2]
	tau  =  P[3]

	arr_t0=arr_t
	arr_t0[0]=0
	kep  = ktrans/ve


	RETURN,  D*KTRANS*(a1*(exp(-kep*(arr_t-tau))-exp(-m1*(arr_t-tau)))/(m1-kep)  + $
					   a2*(exp(-kep*(arr_t-tau))-exp(-m2*(arr_t-tau)))/(m2-kep)) + $

					  D*vp*(a1*(exp(-m1*(arr_t0-tau))) + a2*(exp(-m2*(arr_t0-tau))))

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION MPFIT_CT_TOFTS_CONVOL, X, P

; Ajuste del modelo de Tofts para la función 'MPFITFUN'
; Concentración en tejido (Ct)
;
; Convolución directa con la curva Cp, sin estimación analítica de Cp biexponencial
;
; Ajuste de Ktrans y ve, vp y tau_t
; (vp y tau_t se puede dejar fijo como parametro, o ignorarse poniéndolo a Cero, en la función MPFITFUN)
; El shift temporal (tau) (Orton et al, 2007), también se puede ignorar o dejar fijo

; Al ajustar Ve en lugar de kep (MPFIT_CT_TOFTS_CONVOL_V1), se puede limitar a [0,1] y es más robusto

n_x = N_ELEMENTS(X)
n_times = n_x/2

arr_t = X[0:n_x/2-1]
arr_Cp= X[n_x/2:n_x-1]

IF N_ELEMENTS(arr_t) NE N_ELEMENTS(arr_Cp) THEN RETURN, -1

KTRANS= P[0]
ve   =  P[1]
vp   =  P[2]
tau  =  P[3]

kernel = ktrans[0]*EXP(-(ktrans[0]/ve[0])*(arr_t-tau))

arr_zeros = FLTARR(n_times)
arr_Ct    = CONVOL([arr_zeros, arr_Cp], REVERSE(kernel), NORMALIZE=normalize)

;-----------------------------------------------------------------
p1 = n_times - n_times/2      ; cambio para compatibilidad MATLAB (sumar 1)
p2 = n_times + n_times/2 - 1
IF ((n_times MOD 2) EQ 0) THEN BEGIN
		p1++ & p2++ ; cambio para compatibilidad MATLAB (sumar 1) caso par
	ENDIF; caso par
IF ((p2-p1+1) LT n_times) THEN p2++ ; caso de arrays impares
;-----------------------------------------------------------------

arr_Ct = arr_Ct[p1:p2]

norm_time = (ABS(arr_t[0]-arr_t[n_times-1]))/(n_times-1)
arr_Ct*=norm_time

arr_Ct = arr_Ct + vp[0]*arr_Cp

RETURN, arr_Ct

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION MPFIT_CT_TOFTS_HORSFIELD_V1, X, P

; Estimación según ecuación (4) de Horsfield & Morgan, 2004

n_x = N_ELEMENTS(X)
; n_x es impar...
n_times = n_x/2

; el tiempo está normalizado a [0,1,2...]

Ct0     = X[0]
arr_Ct_1= X[0:n_times-2]
arr_Cp  = X[n_times+1:n_x-1]

IF N_ELEMENTS(arr_Ct_1) NE N_ELEMENTS(arr_Cp) THEN RETURN, -1

ktrans_dt = P[0]   ; ktrans*delta_T
;Eparam    = P[1]   ; EXP(-ktrans*delta_t/ve)

ve        = P[1]
Eparam    = EXP(-ktrans_dt/ve)

arr_Ct =  [Ct0, arr_Ct_1*Eparam + Ktrans_dt*arr_Cp]


RETURN, arr_Ct

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION MPFIT_CT_TOFTS_HORSFIELD_V2, X, P

; Estimación según ecuación (6) de Horsfield & Morgan, 2004

n_x = N_ELEMENTS(X)
; n_x es impar...
n_times = n_x/2

; el tiempo está normalizado a [0,1,2...]

Ct0     = X[0]
arr_Ct_1= X[0:n_times-2]
arr_Cp  = X[n_times+1:n_x-1]
arr_Cp_1= X[n_times:n_x-2]

IF N_ELEMENTS(arr_Ct_1) NE N_ELEMENTS(arr_Cp) THEN RETURN, -1

ktrans_dt = P[0]   ; ktrans*delta_T
;Eparam    = P[1]   ; EXP(-ktrans*delta_t/ve)
ve        = P[1]
Eparam    = EXP(-ktrans_dt/ve)

arr_Ct =  [Ct0, arr_Ct_1*Eparam + Ktrans_dt*(  (arr_Cp*(SQRT(Eparam)-1) + arr_Cp_1*(Eparam - SQRT(Eparam)))/ALOG(Eparam))]

RETURN, arr_Ct

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION MPFIT_CT_TOFTS_HORSFIELD_V3, X, P

; Estimación según ecuación (6) de Horsfield & Morgan, 2004

n_x = N_ELEMENTS(X)
; n_x es impar...
n_times = n_x/2

; el tiempo está normalizado a [0,1,2...]

Ct0     = X[0]
arr_Ct_1= X[0:n_times-2]
arr_Cp  = X[n_times+1:n_x-1]
arr_Cp_1= X[n_times:n_x-2]

IF N_ELEMENTS(arr_Ct_1) NE N_ELEMENTS(arr_Cp) THEN RETURN, -1

ktrans_dt = P[0]   ; ktrans*delta_T
;Eparam    = P[1]   ; EXP(-ktrans*delta_t/ve)

ve        = P[1]
Eparam    = EXP(-ktrans_dt/ve)

arr_Ct =  [Ct0, arr_Ct_1*Eparam + Ktrans_dt*((arr_Cp*(Eparam-ALOG(Eparam)-1) -  arr_Cp_1*(Eparam-Eparam*ALOG(Eparam)-1))/(ALOG(Eparam)^2))]

RETURN, arr_Ct

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION MPFIT_CT_YANKEELOV_CONVOL_V1, X, P

; Ajuste del modelo RR de Yankeelov para la función 'MPFITFUN'
; Concentración en tejido (Ct)
;
; Ajuste de la versión discreta (convolución) de la ecuación (6) en Barnes et al (Pharmaceutics, 2012)
;
; Ajuste de Ktrans y Ve, conocidos Ktrans_RR, Ve_RR (vp_t y vpRR se asumen cero)
;

bg = 2
n_times = (N_ELEMENTS(X)-2)/2

ktrans_RR = X[0]
Ve_RR     = X[1]
arr_time  = X[bg:n_times+bg-1]
arr_Crr   = X[n_times+bg:n_times*2+bg-1]

IF N_ELEMENTS(arr_t) NE N_ELEMENTS(arr_Cp) THEN RETURN, -1

ktrans_t= P[0]
ve_t   =  P[1]

kernel = (ktrans_T/ktrans_RR)*(ktrans_RR/ve_RR - ktrans_t/ve_t)*EXP(-(ktrans_t/ve_t)*arr_time)

arr_zeros = FLTARR(n_times)
arr_Ct    = CONVOL([arr_zeros, arr_Crr], REVERSE(kernel), NORMALIZE=normalize)

;-----------------------------------------------------------------
p1 = n_times - n_times/2      ; cambio para compatibilidad MATLAB (sumar 1)
p2 = n_times + n_times/2 - 1
IF ((n_times MOD 2) EQ 0) THEN BEGIN
		p1++ & p2++ ; cambio para compatibilidad MATLAB (sumar 1) caso par
	ENDIF; caso par
IF ((p2-p1+1) LT n_times) THEN p2++ ; caso de arrays impares
;-----------------------------------------------------------------

arr_Ct = arr_Ct[p1:p2]

norm_time = (ABS(arr_time[0]-arr_time[n_times-1]))/(n_times-1)
arr_Ct*=norm_time

arr_Ct = arr_Ct + (ktrans_t/ktrans_RR)*arr_Crr

RETURN, arr_Ct

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>



FUNCTION MPFIT_RCE_HOFFMANN_V1, X, P

; Función de Hoffmann model adecuado para MPFIT
; (formula en Hoffmann et al, 1995, Tofts, 1997, Huang et al. 2004)
; Notación: (Tofts, 1997)

	n_x = N_ELEMENTS(X)

	arr_t = X[0:n_x-1]

	kep = P[0]
	kel = P[1]
	A   = P[2]

	RETURN,  (A*(kep*(exp(-kep*arr_t)-exp(-kel*arr_t))/(kel-kep))) + 1

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION MPFIT_RCE_HOFFMANN_V2, X, P

; Función de Hoffmann model adecuado para MPFITFUN
; (formula en Hoffmann et al, 1995, Huang et al. 2004)
; Notación: (Tofts, 1997)
;
; Estiamción de RCE
;
; Versión modificada para estimación directa de Akep (en lugar a A)

	n_x = N_ELEMENTS(X)

	arr_t = X[0:n_x-1]

	kep   = P[0]
	kel   = P[1]
	Akep  = P[2]

	RETURN,  (Akep*(exp(-kep*arr_t)-exp(-kel*arr_t))/(kel-kep)) + 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION MPFIT_RCE_HOFFMANN_V3, X, P

; Ajuste del modelo de Hoffmann para la función 'MPFITFUN'
; Contraste relativo (RCE)
;
; Ajuste de kep, kel y A y tiempo de ajuste (shift) tau

; Equivalente al ajuste para ktrans de (Orton et al, 2007)

	n_x = N_ELEMENTS(X)

	arr_t = X[0:n_x-1]

	kep = P[0]
	kel = P[1]
	A   = P[2]
	tau = P[3]

	RETURN,  (A*(kep*(exp(-kep*(arr_t-tau))-exp(-kel*(arr_t-tau)))/(kel-kep))) + 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION MPFIT_RCE_HOFFMANN_V4, X, P

; Ajuste del modelo de Hoffmann para la función 'MPFITFUN'
; Contraste relativo (RCE)
;
; Ajuste de kep, kel y A.

; (formula mas complicada con tiempo de infusión en Hoffmann et al, 1995 (tau_inf)
; Notación: (Tofts, 1997)
;
; Es el modelo de Brix (tau infusion no igual a cero)

	n_x = N_ELEMENTS(X)

	tau_inf = X[0]
	arr_t = X[1:n_x-1]

	kep = P[0]
	kel = P[1]
	A   = P[2]

	rel1 = kep/(kel*(kep-kel))
	rel2 = 1/(kep-kel)

	arr_tp = (arr_t < tau_inf) ; during infusion (t < tau), tp=t, after infusion, tp=tau

	RETURN,  ((A/tau_inf)*(rel1*(exp(kel*arr_tp)-1)*exp(-kel*arr_t)  - $
		                   rel2*(exp(kep*arr_tp)-1)*exp(-kep*arr_t) ) + 1)
END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION MPFIT_RCE_HOFFMANN_V5, X, P

; Ajuste del modelo de Hoffmann para la función 'MPFITFUN'
; Contraste relativo (RCE)
;
; Modificación nueva de modelo de Hoffmann con término adicional "Vp*Cp" pero con cte B
;
; Ajuste de kep, kel y A , B


	n_x = N_ELEMENTS(X)

	arr_t   = X[0:n_x/2-1]
	arr_aif = X[n_x/2:n_x-1]
	; la AIF debe de ser cero al origen, luego subir...

	kep = P[0]
	kel = P[1]
	A   = P[2]
	B   = P[3]

   	RETURN,  (A*(kep*(exp(-kep*arr_t)-exp(-kel*arr_t))/(kel-kep))) + 1 + B*arr_aif

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION MPFIT_RCE_LARSSON_V1, X, P

; Ajuste del modelo de Larsson para la función 'MPFITFUN'
; Contraste relativo (RCE)
;
; Ajuste de Si u kep

n_x = N_ELEMENTS(X)

a1    = X[0]
a2    = X[1]
m1    = X[2]
m2    = X[3]
arr_t   = X[4:n_x-1]

kep = P[0]
Si  = P[1]

RETURN,  1.0 + (Si/(a1+a2)) * ( a1*(EXP(-kep*arr_t) - EXP(-m1*arr_t))/(m1-kep) + $
					      a2*(EXP(-kep*arr_t) - EXP(-m2*arr_t))/(m2-kep) )

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION NormalizedGamma_variate, arr_t, ALPHA=alpha, TAU=tau

;---------------------------------------------------------------------------------
;M. C. Schabel, and D. L. Parker,
;Uncertainty and bias in contrast concentration measurements using spoiled gradient echo pulse sequences,
;Physics in Medicine and Biology, vol. 53, no. 9, pp. 2345-2373, May 7, 2008.
;---------------------------------------------------------------------------------
; (EQUATION A.2)


;arr_result = arr_t^(alpha[0])*EXP(-arr_t/beta[0])/(beta^(alpha+1)*GAMMA(alpha+1))

; Double precission
alphad = DOUBLE(alpha[0])
taud   = DOUBLE(tau[0])

pos_1 = WHERE(arr_t GE 0, ct1)
pos_2 = WHERE(arr_t LT 0, ct2)

arr_result = DBLARR(N_ELEMENTS(arr_t))

IF ct1 GT 0 THEN BEGIN
	arr_t_1    = arr_t[pos_1]

	;arr_result_1 = ((EXP(1)/GAMMA(alpha[0]-1))^(alpha[0]-1))* $
	arr_result_1 = ((EXP(1)/(taud*(alphad-1)))^(alphad-1))* $
			(arr_t_1^(alphad-1)) * $
			(EXP(-arr_t_1/taud))

	arr_result[pos_1] = arr_result_1
ENDIF
IF ct2 GT 0 THEN BEGIN
	arr_result[pos_2] = 0
ENDIF

RETURN, REFORM(arr_result)

END


;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Plot_profile, array_y, array_x, WIN_NUMBER=win_number, TITLE=title, XTITLE=xtitle, YTITLE=ytitle,$
	POINT=point, HR_POINT=hr_point, NO_PLOT_PSYM=no_plot_psym, NO_PLOT_LINE=no_plot_line

LOADCT, 0
DEVICE, DECOMPOSED=1
DEVICE, Get_Screen_Size = screenSize

str_point = ''
IF N_ELEMENTS(point) EQ 2 THEN BEGIN
	 str_point = ' at ROI point: [' + STRTRIM(point[0],2) + ', ' + STRTRIM(point[1],2) + ']
ENDIF

IF N_ELEMENTS(hr_point) EQ 2 THEN BEGIN
	 str_point = ' at MR point: [' + STRTRIM(hr_point[0],2) + ', ' + STRTRIM(hr_point[1],2) + ']
ENDIF

IF N_ELEMENTS(title) THEN str_title=title[0] ELSE str_title=''

IF KEYWORD_SET(no_plot_psym) THEN opt_psym = 0 ELSE opt_psym = 1
IF KEYWORD_SET(no_plot_line) THEN opt_line = 0 ELSE opt_line = 1

WINDOW, win_number, XSIZE=700, YSIZE=500, XPOS=screenSize[0]-720, YPOS=20

yrange = [MIN(array_y)*0.9, MAX(array_y)*1.05]

IF N_ELEMENTS(array_x) EQ 0 THEN BEGIN
	arr_x = INDGEN(N_ELEMENTS(array_y))
	xrange = [0, N_ELEMENTS(arr_x)]

	PLOT,  [0],[0], YRANGE=yrange, XRANGE=xrange, $
		BACKGROUND='FFFFFF'X, COLOR='00'x, XSTYLE=1,YSTYLE=1,$
		XTITLE=xtitle,YTITLE=ytitle, FONT=0, $
		TITLE= str_title + str_point

	IF opt_line THEN $
		OPLOT, arr_x, array_y, COLOR='FF0000'x, THICK=1
	IF opt_psym THEN $
		OPLOT, arr_x, array_y, COLOR='FF0000'x, THICK=2, PSYM=4

ENDIF ELSE BEGIN
	IF N_ELEMENTS(array_x) NE N_ELEMENTS(array_x) THEN RETURN, -1
	xrange = [0, MAX(array_x)]

	PLOT,  [0],[0], YRANGE=yrange, XRANGE=xrange, $
		BACKGROUND='FFFFFF'X, COLOR='00'x, XSTYLE=1,YSTYLE=1,$
		XTITLE=xtitle, YTITLE=ytitle, FONT=0, $
	   	TITLE= str_title + str_point

	IF opt_line THEN $
		OPLOT, array_x, array_y, COLOR='FF0000'x, THICK=1
	IF opt_psym THEN $
		OPLOT, array_x, array_y, COLOR='FF0000'x, THICK=2, PSYM=4

ENDELSE

RETURN, 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Plot_profile_dynMRI, arr_dynMRI, ARR_TIME=arr_time, WIN_NUMBER=win_number, $
	POINT=point, HR_POINT=hr_point, FRAME_INJECTION=frame_injection, METHOD=method, XPOS=xpos, YPOS=ypos


; METHOD=
;   1 - (or not defined) Calculates values at %100 with respect to the MEAN of the baseline (0:frame_injection) values
;   2 - Calculates values at %100 with respect to the first frame
;   3 - Does nothing


IF N_ELEMENTS(method) EQ 0 THEN opt_method=1l ELSE opt_method = FIX(method)

IF N_ELEMENTS(frame_injection) EQ 0 THEN frame_injection = 0l

str_point = ''

CASE opt_method OF
	1 : BEGIN
		S_o = MEAN(arr_dynMRI[0:frame_injection])
		arr_plot=arr_dynMRI/S_o
		arr_plot*=100
		str_title  = 'Signal (% relative to frames before injection)'
		str_ytitle = '%'
	END
	2 : BEGIN
		arr_plot = arr_dynMRI/arr_dynMRI[0]
		arr_plot*=100
		str_title  = 'Signal (% relative to first frame)'
		str_ytitle = '%'
	END
	3 : BEGIN
		arr_plot = arr_dynMRI
		str_title  = 'Signal ' + str_point
		str_ytitle = ''
	END
	ELSE: RETURN, -1
ENDCASE

LOADCT, 0
DEVICE, DECOMPOSED=1
DEVICE, Get_Screen_Size = screenSize

IF N_ELEMENTS(point) EQ 2 THEN BEGIN
	 str_point = ' at ROI point: [' + STRTRIM(point[0],2) + ', ' + STRTRIM(point[1],2) + ']
ENDIF
IF N_ELEMENTS(hr_point) EQ 2 THEN BEGIN
	 str_point = ' at MR point: [' + STRTRIM(hr_point[0],2) + ', ' + STRTRIM(hr_point[1],2) + ']
ENDIF

xsize = 700
ysize = 500
IF N_ELEMENTS(xpos) EQ 0 THEN xpos= screenSize[0]-xsize-20
IF N_ELEMENTS(ypos) EQ 0 THEN Ypos= 20
str_title += str_point


IF N_ELEMENTS(arr_time) EQ 0 THEN BEGIN
	arr_x = INDGEN(N_ELEMENTS(arr_dynMRI))
	str_xtitle = 'Frame number'
ENDIF ELSE BEGIN

	IF N_ELEMENTS(arr_time) NE N_ELEMENTS(arr_dynMRI) THEN RETURN, -1
	arr_x = arr_time
	str_xtitle = 'Time (min)'

ENDELSE

colors = Get_color([56,56,200], /HEX)  ; azul bonito
thicks    = 1     ; cambiar con muchas frames

IF N_ELEMENTS(arr_plot) GT 200 THEN type_data=1l ELSE type_data=3l

;DEVICE, SET_FONT = "Arial*20"
str_set_font = "Arial*20"
;---------------------------------------------------------------------------------
ok = PlotProfile(DATA_X=arr_x, DATA_Y=arr_plot, WIN_NUMBER=win_number, $
		TITLE=str_title, XTITLE=str_xtitle, YTITLE=str_ytitle, XRANGE=xrange, YRANGE=yrange, $
		STR_INFO=str_info, XPOS=xpos, YPOS=ypos, XSIZE=xsize,YSIZE=ysize,SET_FONT=str_set_font, $
		TYPE_DATA=type_data, COLORS=colors, THICKS=thicks)
;---------------------------------------------------------------------------------

RETURN, 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Plot_profile_KineticData, data_kinetics, WIN_NUMBER=win_number, NOTITLE=notitle, $
	NO_PLOT_PSYM=no_plot_psym,   NO_PLOT_LINE=no_plot_line, PLOT_ERROR=plot_error, $
	PLOT_BASELINE=plot_baseline, PLOT_INJECTION=plot_injection,XPOS=xpos, YPOS=ypos, $
	TYPE_DATA=type_data, INJECTION_POINT=injection_point, SCALE_FACTOR=scale_factor


DEVICE, Get_Screen_Size = screenSize
xsize = 700
ysize = 500
IF N_ELEMENTS(xpos) EQ 0 THEN xpos = screenSize[0]-xsize-20
IF N_ELEMENTS(Ypos) EQ 0 THEN ypos = 20

str_point = ''

arr_time  = REFORM(data_kinetics[*,0])
arr_mean  = REFORM(data_kinetics[*,1])
arr_stdev = REFORM(data_kinetics[*,2])

IF KEYWORD_SET(no_plot_psym) THEN type_data = 1
IF KEYWORD_SET(no_plot_line) THEN type_data = 2

n_points = N_ELEMENTS(arr_time)
IF n_points LT 2 THEN RETURN,-1

IF KEYWORD_SET(plot_error) THEN BEGIN
	yrange = [MIN(arr_mean-ABS(arr_stdev/2)), MAX(arr_mean+ABS(arr_stdev/2))]
ENDIF ELSE BEGIN
	yrange = [MIN(arr_mean), MAX(arr_mean)]
ENDELSE
margin = (yrange[1]-yrange[0])/20.0
yrange[0]-=margin
yrange[1]+=margin


xrange = [arr_time[0]-arr_time[1]/2, arr_time[n_points-1]+arr_time[1]/2]


IF KEYWORD_SET(notitle) EQ 0 THEN $
str_title = 'Dynamic MRI (% Relative to mean  before injection) ' + str_point
xtitle    = 'Time (min)'

;type_data = 3
thicks = [1]

colors1 = Get_color([56,56,200], /HEX)  ; azul bonito
colors2 = Get_color([46,120,46], /HEX)  ; azul bonito
color_labels = '000000'X
color_titles = '000000'X

;DEVICE, SET_FONT = "Arial*20"
str_set_font = "Arial*16"
position = [0.05, 0.05, 0.90, 0.95]
;---------------------------------------------------------------------------------
ok = PlotProfile(DATA_X=arr_time, DATA_Y=arr_mean, WIN_NUMBER=win_number, $
		TITLE=str_title, XTITLE=str_xtitle, YTITLE=str_ytitle, XRANGE=xrange, YRANGE=yrange, $
		STR_INFO=str_info, XPOS=xpos, YPOS=ypos, XSIZE=xsize,YSIZE=ysize,POSITION=position,$
		TYPE_DATA=type_data, COLORS=colors1, THICKS=thicks, SET_FONT=str_set_font)
;---------------------------------------------------------------------------------
IF KEYWORD_SET(plot_baseline) THEN BEGIN
	OPLOT, xrange, [100,100], COLOR=color_labels, LINESTYLE=2, THICK=0.5
ENDIF
IF KEYWORD_SET(plot_injection) THEN BEGIN
ok = PlotProfile(DATA_X=arr_time[0:injection_point], DATA_Y=arr_mean[0:injection_point], $
		WIN_NUMBER=win_number, OPLOTT=1, $
		TYPE_DATA=type_data, COLORS=colors2, THICKS=thicks, SET_FONT=str_set_fonf)
ENDIF
IF KEYWORD_SET(plot_error) THEN BEGIN
	ERRPLOT, arr_time, arr_mean-arr_stdev/2.0 , arr_mean+arr_stdev/2, $
		COLOR=colors1, LINESTYLE=0
ENDIF

;pos_title_y1 = [0.03, (!Y.WINDOW[1] - !Y.WINDOW[0])/ 2 + !Y.WINDOW[0]]
;pos_title_y2 = [0.97, (!Y.WINDOW[1] - !Y.WINDOW[0])/ 2 + !Y.WINDOW[0]]
;XYOUTS, pos_title_y1[0],pos_title_y1[1], 'Relative to mean  before injection (%)',$
;	FONT=0, ALIGNMENT=0.5, TEXT_AXES=2, /NORMAL, COLOR=color_titles, ORIENTATION=270
IF N_ELEMENTS(scale_factor) NE 0 THEN BEGIN
	AXIS, YAXIS=1, YRANGE = (!Y.CRANGE)*scale_factor, YSTYLE = 1, $
   		COLOR=color_labels, FONT=0, XTICKFORMAT='(F12.0)', YTICKLAYOUT =1
   	;XYOUTS, pos_title_y2[0],pos_title_y2[1], 'Absolute values (a.u)', FONT=1, ALIGNMENT=0.5,$
   	;	/NORMAL,  COLOR=color_titles, ORIENTATION=90
ENDIF

RETURN, 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Plot_xyouts_matrix, matrix_text, XPOS=xpos, YPOS=ypos, COLOR=color, FONT=font,$
	SEP_COL=sep_col, SEP_ROW=sep_row, SET_FONT=set_font, BOLD=bold

tt_font = 0

IF N_ELEMENTS(set_font) NE 0 THEN BEGIN
	str_set_font = set_font
	opt_bold = KEYWORD_SET(bold)
ENDIF ELSE BEGIN
	opt_bold = 0
ENDELSE

size_matrix = SIZE(matrix_text, /DIMENSIONS)
n_rows  = size_matrix[0]
n_lines = size_matrix[1]

;xsize_win = !D.x_size
;ysize_win = !D.y_size
;char_size_y = !D.y_ch_size
;char_size_x = !D.x_ch_size


IF N_ELEMENTS(sep_row) EQ 0 THEN sep_row = 1.8

width = sep_row[0]* (!D.y_ch_size*1.0/!D.y_size) ; Espaciado y medio

arr_x = FLTARR(n_lines)+ xpos
arr_y = REVERSE(INDGEN(n_lines))*width
arr_y+= (ypos-MAX(arr_y))

IF N_ELEMENTS(sep_col) EQ 0 THEN sep_col = 3.0

sep_additional = FLOAT(sep_col[0])*!D.x_ch_size/!D.x_size

alignment=0
FOR i=0l, n_rows-1 DO BEGIN
	width_max = 0
	FOR j=0l, n_lines-1 DO BEGIN
		IF N_ELEMENTS(set_font) NE 0 THEN BEGIN
			;-------------------------------------------
			IF opt_bold THEN BEGIN
				IF  i EQ 0 OR j EQ 0 THEN BEGIN
					DEVICE, SET_FONT= set_font + '*bold';, TT_FONT=tt_font
				ENDIF ELSE BEGIN
					DEVICE, SET_FONT= set_font
				ENDELSE
			ENDIF ELSE BEGIN
				DEVICE, SET_FONT= set_font;, TT_FONT=tt_font
			ENDELSE
			;-------------------------------------------
		ENDIF

		XYOUTS, arr_x[j], arr_y[j], matrix_text[i,j], /NORMAL, COLOR=color, FONT=font, WIDTH=width, ALIGNMENT=alignment
			width_max = width_max > width
	ENDFOR
	arr_x+=(width_max+sep_additional)
ENDFOR

RETURN, 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION PlotProfile, DATA_X=data_x, DATA_Y=data_y, WIN_NUMBER=win_number, $
	ADD=add, OPLOTT=oplott, COLOR_BACKGROUND=color_background, $
	TITLE=title, XTITLE=xtitle, YTITLE=ytitle, TYPE_DATA=type_data, STR_INFO=str_info,$
	XRANGE=xrange, YRANGE=yrange, SET_FONT=set_font,$
	COLORS=colors, THICKS=thicks, $
	XSIZE=xsize, YSIZE=ysize, YPOS=ypos, XPOS=xpos, _EXTRA=extra

; Type data:
; 1 - Linea continua
; 2 - Puntos sin linea
; 3 - Linea + puntos

IF N_ELEMENTS(win_number) EQ 0 THEN win_number = (!D.WINDOW+1) mod 32

opt_Add = KEYWORD_SET(add) ; continue drawing in the current window

opt_oplott = KEYWORD_SET(oplott) ; only oplot in the same plot

opt_high_quality = 0

IF N_ELEMENTS(xsize) EQ 0 THEN xsize = 800l
IF N_ELEMENTS(ysize) EQ 0 THEN ysize = 600l
IF N_ELEMENTS(xpos) EQ 0 THEN  xpos = 20l
IF N_ELEMENTS(ypos) EQ 0 THEN  ypos = 20l

IF opt_oplott EQ 0 THEN BEGIN

	IF opt_add EQ 0 THEN BEGIN
		LOADCT, 0
		DEVICE, DECOMPOSED=1
		IF N_ELEMENTS(color_background) EQ 0 THEN $
			color_background = Get_color([234,234,234], /HEX)  ; gris claro
		color_Axis = '00'x
		WINDOW, win_number, XSIZE=xsize, YSIZE=ysize, XPOS=xpos, YPOS=ypos
	ENDIF ELSE BEGIN
		IF N_ELEMENTS(color_background) EQ 0 THEN $
			color_background = Get_color([234,234,234], /HEX)  ; gris claro
		color_Axis = '00'x
	ENDELSE
ENDIF

ok = psym_usersim(30, FILL=1)
opt_symbol = 8; user defined with psym_usersim


IF N_ELEMENTS(title)  THEN str_title=title[0]   ELSE str_title=''
IF N_ELEMENTS(xtitle) THEN str_xtitle=xtitle[0] ELSE str_xtitle=''
IF N_ELEMENTS(ytitle) THEN str_ytitle=ytitle[0] ELSE str_ytitle=''

n_points    = N_ELEMENTS(data_x)
size_data_y = SIZE(data_y, /DIMENSIONS)
size_data_x = SIZE(data_x, /DIMENSIONS)

IF N_ELEMENTS(size_data_y) EQ 2 THEN n_plots = size_data_y[1] ELSE n_plots = 1l

IF n_points NE size_data_y[0] THEN RETURN,-1

IF N_ELEMENTS(data_x) EQ 0 THEN BEGIN
	data_xb = INDGEN(n_points)
ENDIF


IF N_ELEMENTS(type_data) EQ 1 THEN type_data = REPLICATE(type_data, n_plots)
IF N_ELEMENTS(type_data) EQ 0 THEN type_data = REPLICATE(1, n_plots)

IF N_ELEMENTS(type_data) NE n_plots THEN RETURN, -2

IF N_ELEMENTS(xrange) NE 2 THEN BEGIN
	diff_x = ABS(data_x[1]-data_x[0])
	xrange = [MIN(data_x)-diff_x , MAX(data_x)+diff_x]
ENDIF
IF N_ELEMENTS(yrange) NE 2 THEN BEGIN
	yrange = [MIN(data_y)*0.9 < MIN(data_y)*1.1 , MAX(data_y)*1.05]
ENDIF
;-----------------------------------------------
IF N_ELEMENTS(thicks) EQ 0 THEN BEGIN
	arr_thicks = REPLICATE(1, n_plots)
ENDIF ELSE BEGIN
	arr_thicks = FIX(thicks)
ENDELSE
;-----------------------------------------------
IF N_ELEMENTS(colors) EQ 0 THEN BEGIN
	arr_colors = LONARR(n_plots)
	arr_colors[0] = Get_color([200,0,0], /HEX)    ; rojo
	IF n_plots GT 1 THEN $
		arr_colors[1] = Get_color([56,56,200], /HEX)  ; azul bonito
	IF n_plots GT 2 THEN $
		arr_colors[2] = Get_color([0,128,0], /HEX)    ; verde oscuro
ENDIF ELSE BEGIN
	arr_colors = colors
ENDELSE
;-----------------------------------------------

;DEVICE, SET_FONT = "Garamond*ITALIC*24"
;DEVICE, SET_FONT = "Helvetica*20"
;DEVICE, SET_FONT = "Arial*20"

IF opt_oplott eq 0 THEN BEGIN

	IF N_ELEMENTS(set_font) NE 0 THEN BEGIN
		DEVICE, SET_FONT=set_font
	ENDIF

	PLOT,  [0],[0], YRANGE=yrange, XRANGE=xrange, CHARSIZE=rel_resol, $
		BACKGROUND=color_background, COLOR=color_axis, XSTYLE=1,YSTYLE=1, NOERASE=opt_add, $
		XTITLE=str_xtitle, YTITLE=str_ytitle, FONT=0,TITLE= str_title, _EXTRA=extra
ENDIF

;xloc = !X.Window[0] - 0.075
;yloc = (0.9 - 0.15) / 2.0 + 0.15
;XYOutS, xloc, yloc, str_ytitle, ALIGNMENT =0.5,ORIENTATION=90.0, TEXT_AXES=3, /NORMAL, COLOR=color_axis, FONT=0

FOR i=0, n_plots-1 DO BEGIN
	CASE type_Data[i] OF
	1 : BEGIN ; linea continua
		OPLOT, data_x, data_y[*,i], COLOR=arr_colors[i], THICK=arr_thicks[i], LINESTYLE=0
	END
	2 : BEGIN
		OPLOT, data_x, data_y[*,i], COLOR=arr_colors[i], THICK=arr_thicks[i], PSYM=opt_symbol
	END
	3 : BEGIN
		OPLOT, data_x, data_y[*,i], COLOR=arr_colors[i], THICK=arr_thicks[i], LINESTYLE=0
		OPLOT, data_x, data_y[*,i], COLOR=arr_colors[i], THICK=arr_thicks[i], PSYM=opt_symbol
	END
	ELSE:
	ENDCASE
ENDFOR

RETURN, 1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION psym_usersim, n_points, FILL=fill, THICK=thick

IF N_ELEMENTS(n_points) EQ 0 THEN n_points = 30

A = FINDGEN(n_points) * (!PI*2/(n_points-1))
; Define the symbol to be a unit circle with 16 points,
; and set the filled flag:
USERSYM, COS(A), SIN(A), FILL=fill, THICK=thick

RETURN ,1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION RCE_Signal_to_tissueConcentration_Barboriak, Signal, S_o=S_o, T10=t10, TR=tr, R1=r1, ANGLE=angle, INJECTION_FRAME=injection_frame,$
	S_equilibrium=S_equilibrium, EXPONENTIAL=exponential

;Barboriak, et al. Journal of Magnetic Resonance Imaging, vol. 27, (6), pp. 1388-1398, 2008.
;Ecuación (2) y (1) para signal_equilibrium
;Originalmente de:
; Li et al., J Magn Reson Imaging 2000;12:347-357.

;Signal = valores absolutos de señal

;t10,   valor T1 inicial
;TR ,   tiempo de repeticion
;Angle: Flip_angle, en radianes
;r1    ; relaxividad (4.39 1/second/mM at 1.5T)
;
;S_o   ; Baseline de la señal (señal inicial). Se pruede proporcionar el valor medio estimado si la señal es ruidosa


; expon_blood = exp(-1*TR*((contrastR1*CbmilliM)+(1/T1_blood)));
; SI_blood = S0_blood*sin(theta)*(1-expon_blood)/(1-(cos(theta)*expon_blood));
; expon_tissue = exp(-1*TR*((contrastR1*Ct)+(1/T1_tissue)));
; SI_tissue = S0_tissue*sin(theta)*(1-expon_tissue)/(1-(cos(theta)*expon_tissue));


IF N_ELEMENTS(injection_frame) EQ 0 THEN injection_frame=0l; bad situation, there is no mean to S_o_in...


IF N_ELEMENTS(S_o) EQ 0 THEN S_o_in= MEAN(Signal[0:injection_frame]) ELSE S_o_in = S_o[0]

exponential = EXP(-tr[0]/t10[0])

IF N_ELEMENTS(S_equilibrium) EQ 0 THEN BEGIN
	S_equilibrium = S_o_in[0]/((1-exponential)*(SIN(angle)/(1-COS(angle)*exponential)))  ; Equation (1) in Barboriak et al, 2008
ENDIF

temp   = (Signal - S_o_in)/(S_equilibrium*SIN(angle[0])) + (1-exponential)/(1-exponential*COS(angle[0]))

arr_R  = -(1/tr[0])*ALOG( (1 - temp)/(1-(temp*COS(angle[0])))) ; Gadolinium concentration

arr_Ct = (arr_R- MEAN(arr_R[0:injection_frame]))/r1[0]  ;


RETURN , arr_ct


END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION  RCE_signalRel_from_Signal, data, FRAME_INJECTION=frame_injection

CASE SIZE(data, /N_DIMENSIONS) OF
	3 : BEGIN	; Imagen [xsize, ysize, nframes] (10,10,nframes)
		type_data = 'IMAGE'
		xsize = (SIZE(data, /DIMENSIONS))[0]
		ysize = (SIZE(data, /DIMENSIONS))[1]
		nframestotal= (SIZE(data_RCE, /DIMENSIONS))[2]
		n_points    = xsize*ysize
		data_RCE    = FLTARR(xsize,ysize,nframestotal)
	END
 	2 : BEGIN; Array [n_points, n_frames] (100,nframes)
 		type_data   = 'ARRAY'
 		n_points    = (SIZE(data, /DIMENSIONS))[0]
 		nframestotal= (SIZE(data, /DIMENSIONS))[1]
 		data_RCE = FLTARR(n_points, nframestotal)
	END
	1 : BEGIN
		type_data   = 'SIGNAL'
		n_points    = 1
		nframestotal= (SIZE(data, /DIMENSIONS))[0]
		data_RCE = FLTARR(nframestotal)
	END
 	ELSE : RETURN, -1
ENDCASE

FOR n=0l, n_points-1 DO BEGIN
	;----------------------------
	CASE type_data OF
	'IMAGE' : BEGIN
		px = n MOD xsize
		py = n / xsize
		arr_data = REFORM(data[px, py, *])
		;S_rel  = arr_data[frame_injection:frame_injection+n_frames-1]/MEAN(arr_data[0:frame_injection])
		S_rel  = arr_data/MEAN(arr_data[0:frame_injection])
		data_RCE[px,py,*] = S_rel

	END
	'ARRAY' : BEGIN
		arr_data = REFORM(data[n, *])
		;S_rel  = arr_data[frame_injection:frame_injection+n_frames-1]/MEAN(arr_data[0:frame_injection])
		S_rel  = arr_data/MEAN(arr_data[0:frame_injection])
		data_RCE[n,*] = S_rel

	END
	'SIGNAL' : BEGIN
		arr_data = REFORM(data)
		;S_rel  = arr_data[frame_injection:frame_injection+n_frames-1]/MEAN(arr_data[0:frame_injection])
		S_rel  = arr_data/MEAN(arr_data[0:frame_injection])
		data_RCE = S_rel
	END
	ENDCASE
ENDFOR

RETURN, Data_RCE

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Read_AIFfunction, str_file, ARR_TIME=arr_time, ARR_SIGNAL=arr_signal, $
	ERROR_ID=error_id, ERROR_STR=error_str


	error_id = 0
	strarr_data = ReadFile_ASCII(str_file, REMOVE_BLANKS=1, ERROR_ID=error_id, MAX_LINES=1e6, MAX_BYTES=1e7)
	IF N_ELEMENTS(arr_data) EQ 1 THEN BEGIN
		error_str = 'Error reading file' & error_id = -6 & RETURN,-1

	ENDIF
	strarr_data = Remove_Comments(strarr_data, STR_COMMENTS='#')

	n_frames = N_ELEMENTS(strarr_data)
	str_tab  = STRING(BYTE(9B))

	arr_AIF   = DBLARR(2,n_frames)

	FOR i=0l, n_frames-1 DO BEGIN
		floats     = STRSPLIT(STRTRIM(strarr_data[i],2), str_tab+' ', /EXTRACT, COUNT=count) ; O BIEN TABULADOR O BIEN ESPACIO
		IF  i EQ 0 THEN BEGIN
			n_columns = N_ELEMENTS(floats) ; number of elements by column must be equal , 2 or 3
			IF (((n_columns NE 2) OR (n_columns NE 3)) EQ 0) THEN BEGIN
				error_str = 'Error: wrong file format (line + ' + STRTRIM(i,2) + '): ' +  strarr_data[i]
				error_id = 5 & RETURN, -1
			ENDIF
		ENDIF
		IF count NE n_columns THEN BEGIN
			error_str = 'Error: wrong file format (line + ' + STRTRIM(i,2) + '): ' +  strarr_data[i]
			error_id = 5 & RETURN, -1
		ENDIF
		arr_aif[*,i]  = DOUBLE(floats[0:1])
	ENDFOR

	arr_signal = REFORM(arr_aif[1,*])
	arr_time   = REFORM(arr_aif[0,*])

	RETURN, arr_aif

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Read_binary_BrukerData, str_file_imag,$
	ROTATION=rotation,    $
	DIMENSIONS=dimensions,$
	TYPE_DATA=type_data,$
	SWAPENDIAN=swapendian

ok = FILE_TEST(str_file_imag)
IF ok NE 1 THEN RETURN, -1

size_bytes = (FILE_INFO(str_file_imag)).size

CASE N_ELEMENTS(dimensions) OF
4 : BEGIN
	xsize   = dimensions[0]*1l
	ysize   = dimensions[1]*1l
	zsize   = dimensions[2]*1l
	nframes = dimensions[3]*1l
END
3 : BEGIN
	xsize   = dimensions[0]*1l
	ysize   = dimensions[1]*1l
	zsize   = dimensions[2]*1l
	nframes = 1
END
ELSE: RETURN, -1
ENDCASE

IF N_ELEMENTS(type_data) EQ 0 THEN type_data = 2l

; de momento solo int y long con signo
CASE type_data OF
2 : n_bytes = 2
3 : n_bytes = 4
ELSE: RETURN,-1
ENDCASE



nt = nframes*xsize*ysize*zsize

IF nt NE size_bytes/n_bytes THEN RETURN, -1

;PRINT, nt
;PRINT, size_bytes
;PRINT, size_bytes*1.0/nt

data = MAKE_ARRAY(xsize, ysize, zsize, nframes, TYPE=type_data)
;------------------------------------------
OPENR, lun, str_file_imag, /GET_LUN, SWAP_ENDIAN=swapendian
READU, lun, data
FREE_LUN, lun

opt_rotate = KEYWORD_SET(rotation)
IF opt_rotate EQ 1 THEN BEGIN
	rotation = 7

	IF nframes GT 1 THEN BEGIN
		FOR n=0l, nframes -1 DO BEGIN
			FOR k=0l, zsize-1 DO BEGIN
				data[*,*,k,n] = ROTATE(data[*,*,k,n],7)
			ENDFOR
		ENDFOR
	ENDIF ELSE BEGIN
		FOR k=0l, zsize-1 DO BEGIN
			data[*,*,k] = ROTATE(data[*,*,k,0],7)
		ENDFOR
	ENDELSE
ENDIF

IF nframes EQ 1 THEN data = REFORM(data)

RETURN, data

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Read_bruker_T1map, str_file_imag, ROTATION=rotation, DIMENSIONS=dimensions, $
	DATA_SLOPE=data_slope, DATA_OFFSET=data_offset, ERROR_STR=error_str, ST_DATA=st_data


nimages   = 5l ; por defecto
idx_t1map = 2l ; la tercera imagen es el mapa T1

error_str = ''
opt_rotate = KEYWORD_SET(rotation)

size_file = (FILE_INFO(str_file_imag)).size

xsize   = dimensions[0]*1l
ysize   = dimensions[1]
nslices = dimensions[2]

IF N_ELEMENTS(data_offset) EQ 0 THEN data_offset= FLTARR(nslices)
IF N_ELEMENTS(data_slope) EQ 0 THEN  data_slopet= FLTARR(nslices)+1

IF N_ELEMENTS(data_slope)  NE nslices THEN RETURN,-1
IF N_ELEMENTS(data_offset) NE nslices THEN RETURN,-2

size_must_be = xsize*ysize*nslices*nimages*4

IF size_must_be NE size_file THEN BEGIN
	error_Str = 'File size does not match with dimensions and 5 images by frame'
	RETURN, -1
ENDIF

data = READ_RAW(FILE=str_file_imag, DIMENSIONS=[xsize, ysize, nimages, nslices], TYPE='LONG')

im_map_t1 = FLTARR(xsize, ysize, nslices)

FOR i=0l, nslices-1 DO BEGIN
	im_map_T1[*,*,i] = data[*,*,idx_t1map,i]*data_slope[i] + data_offset[i]
	IF opt_rotate THEN $
		im_map_T1[*,*,i] = ROTATE(im_map_T1[*,*,i],7)
ENDFOR

RETURN, im_map_T1

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Read_BrukerBiospin_image, str_file, ST_DATA=st_data, $
	ERROR_ID=error_id, ERROR_STR=error_str, ONLY_HEADER=only_header, NOSLOPE=noslope,$
	FILE_METHODS=file_methods,$
	NFRAMES=nframes, NORECONS=norecons


; Four files that MUST exist
str_file_2dseq  =  '2dseq'
str_file_visu   =  'visu_pars'
str_file_method =  'method'
str_file_reco   =  'reco'

error_id   =-1
error_str = ''

opt_reco_file = KEYWORD_SET(norecons) EQ 0
;----------------------------------------------------------------------------------------------
IF N_ELEMENTS(str_file) NE 1 THEN RETURN, -1
IF FILE_TEST(str_file, /DIRECTORY) EQ 1 THEN BEGIN ; el el directorio

	str_fileseq = FILE_SEARCH(get_complete_path(str_file) + str_file_2dseq, COUNT=count, /TEST_REGULAR)
	IF count NE 1 THEN BEGIN
		error_str ='raw Data file not found (' + str_file_2dseq + ')'
		error_id=-1 & RETURN,-1
	ENDIF
	IF FILE_TEST(str_fileseq[0], /REGULAR) NE 1 THEN BEGIN
		error_str ='raw Data file is not a valid file (' +  str_file_2dseq + ')'
		error_id=-2 & RETURN,-1
	ENDIF
ENDIF ELSE BEGIN
	str_fileseq = str_file
ENDELSE
str_fileseq = str_fileseq[0]
;----------------------------------------------------------------------------------------------
path_1 = get_path(str_fileseq)
str_visu_file = FILE_SEARCH(path_1 + str_file_visu, COUNT=count, /TEST_REGULAR)
IF count EQ 0 THEN BEGIN
	str_visu_file = FILE_SEARCH(path_1 + str_file_visu + '.txt', COUNT=count, /TEST_REGULAR)
ENDIF
IF count NE 1 THEN BEGIN
	error_str = str_file_visu + ' header file not found' & error_id=-3 & RETURN,-1
ENDIF
;----------------------------------------------------------------------------------------------

IF opt_reco_file THEN BEGIN
	str_reco_file = FILE_SEARCH(path_1 + str_file_reco, COUNT=count, /TEST_REGULAR)
	IF count EQ 0 THEN BEGIN
		str_reco_file = FILE_SEARCH(path_1 + str_file_reco + '.txt', COUNT=count, /TEST_REGULAR)
	ENDIF
	IF count NE 1 THEN BEGIN
		error_str = str_file_reco + ' header file not found' & error_id=-3 & RETURN,-1
	ENDIF
ENDIF
;----------------------------------------------------------------------------------------------

IF N_ELEMENTS(file_methods) EQ 0 THEN BEGIN
	path_2 = FILE_DIRNAME(FILE_DIRNAME(path_1), /MARK_DIRECTORY)
	str_methods_file = FILE_SEARCH(path_2 + str_file_method, COUNT=count, /TEST_REGULAR)
	IF count EQ 0 THEN BEGIN
		str_methods_file = FILE_SEARCH(path_2 + str_file_method + '.txt', COUNT=count, /TEST_REGULAR)
	ENDIF
	IF count EQ 0 THEN BEGIN ; en el mismo directorio
		str_methods_file = FILE_SEARCH(path_1 + str_file_method + '.txt', COUNT=count, /TEST_REGULAR)
	ENDIF
	IF count EQ 0 THEN BEGIN ; en el mismo directorio
		str_methods_file = FILE_SEARCH(path_1 + str_file_method, COUNT=count, /TEST_REGULAR)
	ENDIF
	IF count NE 1 THEN BEGIN
		error_str ='method header file not found' & error_id=-4 & RETURN,-1
	ENDIF
ENDIF ELSE BEGIN
	str_methods_file = file_methods[0]
	IF FILE_TEST(str_methods_file) EQ 0 THEN BEGIN
		error_str ='method header file not found' & error_id=-4 & RETURN,-1
	ENDIF

ENDELSE

;----------------------------------------------------------------------------------------------
str_visu_file = str_visu_file[0]
str_methods_file = str_methods_file[0]

IF opt_reco_file THEN str_reco_file = str_reco_file[0]

PRINT, 'Visualization header file :  ' + str_visu_file
PRINT, 'Methods header file:         ' + str_methods_file
PRINT, 'Raw data file:               ' + str_fileseq

IF opt_reco_file THEN BEGIN
	PRINT, 'Reconstructiont header file: ' + str_reco_file
ENDIF

error_id = 0

st_methods = Read_BrukerBiospin_MethodHeader(str_methods_file, ERROR_ID=error_id, ERROR_STR=error_STR)
IF error_id NE 0 THEN BEGIN
	RETURN, -1
ENDIF

st_visu = Read_BrukerBiospin_VisuHeader(str_visu_file, ERROR_ID=error_id, ERROR_STR=error_STR, INFO=info_visu)
IF error_id NE 0 THEN RETURN, -1

IF opt_reco_file THEN BEGIN
	st_reco = Read_BrukerBiospin_RecoHeader(str_reco_file, ERROR_ID=error_id, ERROR_STR=error_STR)
	IF error_id NE 0 THEN RETURN, -1
ENDIF

;--------------------------------------------------------------------------------------------------

IF opt_reco_file THEN BEGIN
	st_data = CREATE_STRUCT(st_methods, st_visu, st_reco, {file:str_fileseq[0], info:info_visu})
ENDIF ELSE BEGIN
	st_data = CREATE_STRUCT(st_methods, st_visu, {file:str_fileseq[0], info:info_visu})
ENDELSE


type_data = st_data.visu_type_data


IF N_ELEMENTS(nframes) EQ 1 THEN BEGIN
	st_data.nframes = nframes ; don't use the value readed in methods
ENDIF


IF KEYWORD_SET(only_header) THEN BEGIN
	RETURN ,1
ENDIF ELSE BEGIN
	xsize = st_data.nfrec
	ysize = st_data.nphases
	data = Read_Binary_BrukerData(str_fileseq, ROTATION=1,$
    	DIMENSIONS=[xsize, ysize, st_data.nslices, st_data.nframes], TYPE_DATA=type_data)
    IF SIZE(data, /N_DIMENSIONS) LE 2 THEN BEGIN
    	error_id = -7
    	error_str = 'Error reading data'
       	RETURN,-1
	ENDIF
	;----------------------------------------------------------------------------------------------
	; multiplica por el slope

	IF KEYWORD_SET(noslope) EQ 0 THEN BEGIN

		IF N_ELEMENTS(st_data.visu_dataslope) EQ 1 THEN BEGIN
			data = FLOAT(data)*st_data.visu_dataslope ; before load image...importante
		ENDIF ELSE BEGIN
			data = FLOAT(data)
			IF st_data.nslices EQ 1 THEN BEGIN
				FOR j=0l, st_data.nframes-1 DO BEGIN
					data[*,*,0,j] = data[*,*,0,j]*st_data.visu_dataslope[j]

				ENDFOR
			ENDIF ELSE BEGIN
				n=0l
				FOR j=0l, st_data.nframes-1 DO BEGIN
					FOR i=0l, st_data.nslices-1 DO BEGIN
						data[*,*,i,j] = data[*,*,i,j]*st_data.visu_dataslope[n]
						n++
					ENDFOR
				ENDFOR
			ENDELSE
		ENDELSE
	ENDIF ELSE BEGIN
		data = FLOAT(data)
	ENDELSE

	RETURN, data
ENDELSE

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Read_BrukerBiospin_MethodHeader, str_file,  $
	ERROR_ID=error_id, ERROR_STR=error_str, T1MAP=t1map


error_id   =-1
error_str = ''

IF KEYWORD_SET(t1map) THEN opt_t1map = 1 ELSE opt_t1map = 0

str_lines = ReadFile_ASCII(str_file, REMOVE_BLANKS=1, ERROR_ID=error_id, $
    MAX_LINES=max_lines, MAX_BYTES=1e5)

IF error_id NE 0 THEN BEGIN
	error_str ='Error reading file' & error_id=-1 & RETURN,-1
	RETURN, -1
ENDIF

n_lines = N_ELEMENTS(str_lines)

n_keys = 5l

str_codes = ['##$PVM_NRepetitions','##$PVM_RepetitionTime', '##$PVM_NAverages',$
	'##$PVM_SPackArrNSlices', '##$PVM_Matrix']

arr_ln = INTARR(n_keys)

FOR j=0l, n_keys-1 DO BEGIN
	pos_strpos = STRPOS(str_lines, str_codes[J])
	pos_exist  = WHERE(pos_strpos NE -1, ct)
	IF ct EQ 1 THEN BEGIN
		arr_ln[j] = pos_Exist[0]
	ENDIF ELSE BEGIN
		IF opt_t1map EQ 1 AND j EQ 1 THEN BEGIN ; special case

		ENDIF ELSE BEGIN
			error_Str = str_codes[j] + ' key not found or duplicated' & error_id = -6 & RETURN, -1
		ENDELSE
	ENDELSE
ENDFOR

FOR j=0l, n_keys-1 DO BEGIN

	IF opt_t1map EQ 1 AND j EQ 1 THEN BEGIN ; special case
		RepetitionTime = 0.0;

	ENDIF ELSE BEGIN
		;----------------------------------------------------------------------------------
		number  = Brukerline_getNumberValue(str_lines[arr_ln[j]], NEXT_LINE=next_line, ERROR_ID=error_id)
		IF error_id NE 0 THEN BEGIN
			error_Str = str_codes[j] + ' wrong format (1)' & error_id = -6 & RETURN, -1
		ENDIF
		IF next_line EQ 1 THEN BEGIN
			number = Brukerline_getValues(str_lines[arr_ln[j]+1], FIX(number),ERROR_ID=error_id)
			IF error_id NE 0 THEN BEGIN
				error_Str = str_codes[j] + ' wrong format (2)' & error_id = -6 & RETURN, -1
			ENDIF
		ENDIF
		;----------------------------------------------------------------------------------
		CASE j OF
		0 : BEGIN
			NRepetitions = LONG(number)
		END
		1 : BEGIN
			RepetitionTime = FLOAT(number)
		END
		2 : BEGIN
			NAverages = LONG(number)
		END
		3: BEGIN
			SPackArrNSlices = LONG(number)
		END
		4: BEGIN
			Matrix = LONG(number)
			IF N_ELEMENTS(matrix) NE 2 THEN RETURN, -1
		END
		ENDCASE
	ENDELSE
ENDFOR
;---------------------------------------------------------------
st_data = {$
	nframes  : NRepetitions[0],$
	RepetitionTime : RepetitionTime[0]/1000.0,$ // en segundos
	NAverages : NAverages[0],$
	nslices :  SPackArrNSlices[0],$
	nphases :  Matrix[1],$
	nfrec   :  matrix[0],$
	frametime : 0d}

st_data.frametime = st_data.RepetitionTime*st_data.NAverages*st_data.nphases


error_id = 0
RETURN, st_data

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Read_BrukerBiospin_RecoHeader, str_file, ERROR_ID=error_id, ERROR_STR=error_STR, INFO=info

error_id   =-1
error_str = ''

; Lee el visuCoreDataslope de la cabecera Bruker

str_array = ReadFile_ASCII(str_file, REMOVE_BLANKS=1, ERROR_ID=error_id, $
    MAX_LINES=1e5, MAX_BYTES=1e6)

IF error_id NE 0 THEN RETURN, -1

pos_info = STRPOS(str_array, '$$')
pos_ss = WHERE(pos_info EQ 0,ct_pos1)
IF ct_pos1 GE 2 THEN BEGIN
	date = STRMID(str_array[pos_ss[0]],3, STRLEN(str_array[pos_ss[0]])-3)
	path = get_path(STRMID(str_array[pos_ss[1]],3, STRLEN(str_array[pos_ss[1]])-3))
	info = [[date], [path]]
ENDIF ELSE BEGIN
    info = ''
ENDELSE

;--------------------------------------------------------------
; los arrays de slopes, datamin, datamax...

reco_minima = GetBrukerHeader_Array(str_array, '##$RECO_minima', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-2
reco_maxima = GetBrukerHeader_Array(str_array, '##$RECO_maxima', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-3
reco_offset = GetBrukerHeader_Array(str_array, '##$RECO_map_offset', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-4
reco_slope  = GetBrukerHeader_Array(str_array,   '##$RECO_map_slope', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-5

;--------------------------------------------------------------
; el  tipo de datos

;--------------------------------------------------------------
; el  tipo de datos
key_Value = GetBrukerHeader_key(str_array, '##$RECO_wordtype', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-6
CASE STRUPCASE(key_Value) OF
'_32BIT_SGN_INT': type_data = 3 ; 'LONG'
'_16BIT_SGN_INT': type_data = 2 ; 'INT'
ELSE: type_data = 4
ENDCASE
;--------------------------------------------------------------

st_data = {reco_dataslope:reco_slope,$
	reco_datamax:reco_maxima, $
	reco_datamin:reco_minima, $
	reco_dataoffs:reco_offset,$
	reco_type_data:type_data  $
}

error_id = 0

RETURN, st_data

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Read_BrukerBiospin_VisuHeader, str_file, ERROR_ID=error_id, ERROR_STR=error_STR, INFO=info

error_id   =-1
error_str = ''

; Lee el visuCoreDataslope de la cabecera Bruker

str_array = ReadFile_ASCII(str_file, REMOVE_BLANKS=1, ERROR_ID=error_id, $
    MAX_LINES=1e5, MAX_BYTES=1e6)

IF error_id NE 0 THEN RETURN, -1

pos_info = STRPOS(str_array, '$$')
pos_ss = WHERE(pos_info EQ 0,ct_pos1)
IF ct_pos1 GE 2 THEN BEGIN
	date = STRMID(str_array[pos_ss[0]],3, STRLEN(str_array[pos_ss[0]])-3)
	path = get_path(STRMID(str_array[pos_ss[1]],3, STRLEN(str_array[pos_ss[1]])-3))
	info = [[date], [path]]
ENDIF ELSE BEGIN
    info = ''
ENDELSE

;--------------------------------------------------------------
; los arrays de slopes, datamin, datamax...
dataslope = GetBrukerHeader_Array(str_array, '##$VisuCoreDataSlope', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-2
datamin   = GetBrukerHeader_Array(str_array, '##$VisuCoreDataMin', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-3
datamax   = GetBrukerHeader_Array(str_array, '##$VisuCoreDataMax', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-4
dataoffs  = GetBrukerHeader_Array(str_array, '##$VisuCoreDataOffs', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-5
;--------------------------------------------------------------
; el  tipo de datos
key_Value = GetBrukerHeader_key(str_array, '##$VisuCoreWordType', ERROR_ID=error_id2)
IF error_id2 NE 0 THEN RETURN,-6
CASE STRUPCASE(key_Value) OF
'_32BIT_SGN_INT': type_data = 3 ; 'LONG'
'_16BIT_SGN_INT': type_data = 2 ; 'INT'
ELSE: type_data = 4
ENDCASE
;--------------------------------------------------------------
error_id = 0
st_data = {visu_dataslope:dataslope,$
	visu_datamin:datamin,$
	visu_datamax:datamax,$
	visu_dataoffs:dataoffs,$
	visu_type_data:type_data $
}

error_id = 0

RETURN, st_data

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>


FUNCTION Read_CpModelConfig, str_file_config, MODEL=model, STRUCT=struct, ERROR_ID=error_id, ERROR_STR=error_str

error_id=0
error_str = ''

array_in  = ReadFile_ASCII(str_file_config, REMOVE_BLANKS=1, ERROR_ID=error_id, MAX_LINES=1e4, MAX_BYTES=1e5)
IF error_id NE 0 THEN BEGIN
	error_str = 'Error reading config file' &	RETURN, -1
ENDIF

array_in = Remove_Comments(array_in, STR_COMMENTS='#')

n_lines_min = 4

n_lines = N_ELEMENTS(array_in)
strarr_lines = array_in
IF n_lines LT n_lines_min THEN BEGIN & error_id = 4 & error_str = 'Error: not valid config file' & RETURN, -1 & ENDIF

;------------------------------------------------------
strarr_values = ReadStrarr_obtainvalues(strarr_lines[0:n_lines-1], STR_SEPARATOR='=', ERROR_ID=err2, ERROR_STR=error_str)
IF err2 NE 0 THEN BEGIN & error_id =7l &  RETURN, -1 & ENDIF
strarr_values[0,*] = STRLOWCASE(strarr_values[0,*])
;------------------------------------------------------
;
string_model = 'model'

pos_model = 0 ; must be the first one
pos = WHERE(STRMATCH(strarr_values[0,pos_model], string_model) EQ 1, ct)
IF ct NE 1 THEN BEGIN
	error_id  = -1 & error_str = 'Keyword <model> not found' & RETURN, -1
ENDIF
var_model = STRLOWCASE((STRTRIM( (strarr_values[1,pos_model])[0],2))[0])

strarr_models  = ['biexponential', 'orton', 'schabel-simplified', 'schabel-standard', 'parker']

IF MAX(var_model EQ strarr_models) NE 1 THEN BEGIN
	error_id  = -1 & error_str = 'Model unknown or not supported:'  + var_model +  ' ?? '& RETURN, -1
ENDIF

CASE var_model OF
	strarr_models[0] :	strarr_Data = ['m1','m2','a1','a2'] ; BIEXPONENTIAL
	strarr_models[1] :	strarr_Data = ['ab','ag','ub','ug'] ; ORTON
	strarr_models[2] :	strarr_Data = ['a1', 'a2', 'alpha', 'tau1', 'tau2', 'delay'] ; schabel-simplified
	strarr_models[3] : 	strarr_Data = ['a1', 'a2', 'a3', 'a4', 'alpha', 'tau1', 'tau2', 'delay', 'd1','d2','d3']  ; schabel-standard
	strarr_models[4] : 	strarr_Data = ['alpha', 'beta', 's', 'tau', 't1', 't2', 'sigma1', 'sigma2', 'a1','a2'] ; parker
ENDCASE
strarr_Data = STRLOWCASE(strarr_data)
;---------------------------------------------
n_data  = N_ELEMENTS(strarr_Data)
;err=0
;ascii_min = BYTE('0')
;ascii_max = BYTE('9')
;
strarr_values = strarr_values[*,1:n_lines-1]
;
n_campos=0l
FOR i=0l, n_data-1 DO BEGIN
	pos = WHERE(STRMATCH(strarr_values[0,*], strarr_data[i]) EQ 1, ct)
	IF ct NE 1 THEN BEGIN
		error_id = -1 &	error_str = 'Keyword not found: ' + strarr_data[i]
		BREAK
	ENDIF
	var_temp = (STRTRIM((strarr_values[1,pos])[0],2))[0]

	CASE strarr_data[i] OF
	 	'm1' : m1 = (FLOAT(var_temp))[0]
	 	'm2' : m2 =(FLOAT(var_temp))[0]
	 	'a1' : a1 =(FLOAT(var_temp))[0]
	 	'a2' : a2 =(FLOAT(var_temp))[0]
	 	'ab' : ab =(FLOAT(var_temp))[0]
	 	'ag' : ag =(FLOAT(var_temp))[0]
	 	'ub' : ub =(FLOAT(var_temp))[0]
	 	'ug' : ug =(FLOAT(var_temp))[0]
	 	'alpha' : alpha = (FLOAT(var_temp))[0]
	 	'tau1' : tau1   = (FLOAT(var_temp))[0]
	 	'tau2' : tau2   = (FLOAT(var_temp))[0]
	 	'delay' : delay = (FLOAT(var_temp))[0]
	 	'a3' : a3 = (FLOAT(var_temp))[0]
	 	'a4' : a4 = (FLOAT(var_temp))[0]
	 	'd1' : d1 = (FLOAT(var_temp))[0]
	 	'd2' : d2 = (FLOAT(var_temp))[0]
	 	'd3' : d3 = (FLOAT(var_temp))[0]
	 	'beta' : beta = (FLOAT(var_temp))[0]
	 	's' : s = (FLOAT(var_temp))[0]
	 	'tau' : tau = (FLOAT(var_temp))[0]
	 	't1' : t1 = (FLOAT(var_temp))[0]
	 	't2' : t2 = (FLOAT(var_temp))[0]
	 	'sigma1' : sigma1 = (FLOAT(var_temp))[0]
	 	'sigma2' : sigma2 = (FLOAT(var_temp))[0]
	ELSE : BEGIN
		error_id = 8l
		error_str = 'not valid keyword : ' + strarr_data[i]
		RETURN,-1
	END
	ENDCASE
    n_campos++
ENDFOR
IF n_campos NE n_data THEN BEGIN & error_id = 9l & RETURN,-1 & ENDIF
;---------------------------------------------------------------------------------

CASE var_model OF
strarr_models[0] : BEGIN
	struct = {model:strarr_models[0], m1:m1,m2:m2,a1:a1,a2:a2}
END
strarr_models[1] : BEGIN
	struct = {model:strarr_models[1], ab:ab,ag:ag,ub:ub,ug:ug}
END
strarr_models[2] : BEGIN
	struct = {model:strarr_models[2], a1:a1, a2:a2, alpha:alpha, tau1:tau1, tau2:tau2, delay:delay}
END
strarr_models[3] : BEGIN
	struct = {model:strarr_models[3], a1:a1, a2:a2, a3:a3, a4:a4, alpha:alpha, tau1:tau1, tau2:tau2, delay:delay, d1:d1,d2:d2,d3:d3}
END
strarr_models[4] : BEGIN
	struct = {model:strarr_models[4], alpha:alpha, beta:beta, s:s, tau:tau, t1:t1,t2:t2, sigma1:sigma1, sigma2:sigma2, a1:a1, a2:a2}
END
ENDCASE

model = var_model
RETURN, struct

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Read_DynamicDicom, PATH=path, ERROR_STR=error_Str, ERROR_ID=error_id

error_id = 1l

arr_files = FILE_SEARCH(Get_complete_path(path) + '*.dcm')
n_files   = N_ELEMENTS(arr_files)

IF n_files eq 1 THEN BEGIN
	IF arr_files[0] EQ '' THEN BEGIN
		ERROR_STR='No DCM files in the path'
		RETURN,-1
	ENDIF
ENDIF
arr_files = arr_files[SORT(arr_files)]
opt_rotate = 7
n_slices = 1
FOR i=0,n_files-1 DO BEGIN

	filename_dicom = arr_files[i]
	IF FILE_TEST(filename_dicom, /REGULAR) NE 1 THEN BEGIN
		error_str = 'Error readinf file' + get_name(filename_dicom)
		RETURN,-2
	ENDIF
	ok = QUERY_DICOM(filename_dicom, Info, IMAGE_INDEX=image_index, DICOMEX=dicomex)

	IF i EQ 0 THEN BEGIN
		xsize = info.dimensions[0]
		ysize = info.dimensions[1]
		im_dym = MAKE_ARRAY(TYPE=info.pixel_type, xsize, ysize, n_slices, n_files)
	ENDIF
	if i MOD 40 EQ 0 THEN PRINT, i, n_files, ok
	IF ok EQ 1 THEN BEGIN
		;HELP, info, /STRUCTURE
		image = READ_DICOM(filename_dicom, IMAGE_INDEX=0)
		im_dym[*,*,0,i] = ROTATE(image, opt_rotate)
		;oDicom = OBJ_NEW('IDLffDICOM' , Filename_dicom, /VERBOSE)
	ENDIF ELSE BEGIN
		PRINT, 'Error reading DICOM image ' + filename_dicom
		RETURN,-3
	ENDELSE
ENDFOR

error_str = ''
error_id=0l
RETURN, im_dym

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
;+
;
;NAME:
;   Read_Raw
;
; -- PROPOSITO
;
;   Lee datos en formato binario (RAW)
;
; -- EJEMPLO DE LLAMADA EN LINEA DE COMANDOS IDL:
;
;
; datos = Read_raw($
;    FILE=file, $
;    OFFSET=offset,
;    DIMENSIONS=dimensions,
;    SWAP=swap,
;    TYPE=type, $
;    REVERSE_XYZ=reverse_xyz,
;    ERROR_ID = error_id,
;    ERROR_STR = error_str)
;
; -- ENTRADAS:
;
;FILE
;   Archivo de lectura de los datos. tipo STRING
;
;OFFSET
;   offset al comienzo del fichero, en número de bytes
;
;DIMENSIONS
;   array de valores con las dimensiones de los datos
;
;SWAP
;   Opción para leer en little endian
;
;TYPE
;   variable de tipo STRING que especifica el tipo de datos
;   Admite los tipos IDL:
;
;	'BYTE'     byte (8 bits)
;	'INT':     Entero corto (2 bytes)
;	'INTEGER': Entero corto (2 bytes)
;	'LONG':    Entero largo (4 bytes)
;	'FLOAT':   Punto flotante de 4 bytes
;	'FLOATING':Punto flotante de 4 bytes
;	'DOUBLE':  Punto flotante de 8 bytes
;	'UINT':    Entero corto sin signo (2 bytes)
;	'ULONG':   Entero largo sin signo (4 bytes)
;	'LONG64':  Entero extralargo (8 bytes)
;	'ULONG64': Entero extralargo sin signo (8 bytes)
;
;-


FUNCTION Read_raw, FILE=file, OFFSET=offset, DIMENSIONS=dimensions, SWAP=swap, TYPE=type, $
    REVERSE_XYZ=reverse_xyz, ERROR_ID = error_id, ERROR_STR = error_str


IF N_ELEMENTS(type) EQ 0 THEN type2 = FIX(4) ELSE type2 = type[0]

IF SIZE(type2, /TYPE) EQ 7 THEN BEGIN ; ?type definido como string
	CASE STRUPCASE(STRTRIM(type2,2)) OF
	'BYTE': 	imag_type = 1
	'INT':  	imag_type = 2
	'INTEGER':  imag_type = 2
	'LONG':  	imag_type = 3
	'FLOAT':  	imag_type = 4
	'FLOATING': imag_type = 4
	'DOUBLE':  	imag_type = 5
	'UINT':  	imag_type = 12
	'UINTEGER': imag_type = 12
	'ULONG':  	imag_type = 13
	'LONG64':  	imag_type = 14
	'ULONG64':  imag_type = 15
	ELSE : RETURN, -1
	ENDCASE
ENDIF ELSE BEGIN
	imag_type=FIX(type2)
ENDELSE

IF N_ELEMENTS(reverse_xyz) NE 3 THEN rev_xyz =[0,0,0b] ELSE rev_xyz=reverse_xyz[0:2]

error_id = 1
error_str = ''

IF N_ELEMENTS(file) NE 1 THEN RETURN, -1
IF SIZE(file, /TYPE) NE 7 THEN RETURN, -1
IF N_ELEMENTS(dimensions) LE 0 THEN RETURN, -1
IF N_ELEMENTS(dimensions) GT 4 THEN RETURN, -1
IF N_ELEMENTS(dimensions) LT 4 THEN n_imags = 1L ELSE n_imags = (LONG(dimensions[3]) > 1)
IF N_ELEMENTS(dimensions) LT 3 THEN zsize = 1L ELSE zsize = (LONG(dimensions[2]) > 1)
xsize = (LONG(dimensions[0]) > 1)
IF N_ELEMENTS(dimensions) GT 1 THEN $
	ysize = (LONG(dimensions[1]) > 1) ELSE ysize = 1

IF N_ELEMENTS(offset) EQ 0 THEN bytes_offset = 0L ELSE bytes_offset = LONG(offset)

IF FILE_TEST(file[0]) NE 1 THEN BEGIN
	error_str = 'File can not be readed' &   error_id = 1
	RETURN, -1
END

arrayAux  = MAKE_ARRAY(xsize, ysize, TYPE=imag_type)
image = MAKE_ARRAY(xsize, ysize, zsize, n_imags, TYPE=imag_type)
image = REFORM(image, xsize, ysize, zsize, n_imags)

; Calculo cual es el tamaño en bytes para asegurarme de que el fichero es del tamaño adecuado
CASE imag_type of
1: size_in_bytes=xsize*ysize*zsize*n_imags  ; BYTE
2: size_in_bytes=2*xsize*ysize*zsize*n_imags;INT
3: size_in_bytes=4*xsize*ysize*zsize*n_imags;LINT
4: size_in_bytes=4*xsize*ysize*zsize*n_imags ;Floating point
5: size_in_bytes=8*xsize*ysize*zsize*n_imags ;Double-precision floating
12: size_in_bytes=2*xsize*ysize*zsize*n_imags;Unsigned Integer
13: size_in_bytes=4*xsize*ysize*zsize*n_imags;Unsigned Longword Integer
14: size_in_bytes=8*xsize*ysize*zsize*n_imags;64-bit Integer
15: size_in_bytes=8*xsize*ysize*zsize*n_imags; Unsigned 64-bit Integer
ENDCASE

OPENR, unit, file[0], /GET_LUN, ERROR=err
IF err NE 0 THEN BEGIN
    IF N_ELEMENTS(unit) GT 0 THEN BEGIN & FREE_LUN, unit & CLOSE, unit & ENDIF
    error_id = 1 & RETURN, -1
ENDIF
IF (FSTAT(unit)).SIZE NE (bytes_offset+size_in_bytes) THEN BEGIN
    error_str = 'File size does not match image dimensions' &   error_id = 1
    FREE_LUN, unit & CLOSE, unit
    RETURN,  -1
ENDIF
POINT_LUN,  unit,  bytes_offset

FOR ni =0, n_imags-1 DO BEGIN
    FOR nz=0, zsize-1 DO BEGIN
       READU, unit, arrayAux
       image[*, *, nz, ni] = arrayAux
    ENDFOR
ENDFOR
FREE_LUN, unit & CLOSE, unit

IF KEYWORD_SET(swap) THEN image=SWAP_ENDIAN(image)

IF n_imags EQ 1 THEN BEGIN
	image = REFORM(image, xsize, ysize, zsize)
	IF  rev_xyz[0] THEN image=REVERSE(image,  1,  /Overwrite)
	IF  rev_xyz[1] THEN image=REVERSE(image,  2,  /Overwrite)
	IF zsize GT 1 THEN BEGIN
		IF rev_xyz[2] THEN image=REVERSE(image,  3,  /Overwrite)
	ENDIF ELSE BEGIN
		image = REFORM(image, xsize, ysize)
	ENDELSE
ENDIF ELSE BEGIN
	IF MAX(rev_xyz) GT 0 THEN BEGIN
		FOR ni=0l, n_imags-1 DO BEGIN
			IF  rev_xyz[0] THEN image[*,*,*,ni] = REVERSE(image[*,*,*,ni],  1,  /Overwrite)
			IF  rev_xyz[1] THEN image[*,*,*,ni] = REVERSE(image[*,*,*,ni],  2,  /Overwrite)
			IF zsize GT 1 THEN BEGIN
				IF rev_xyz[2] THEN  image[*,*,*,ni]=REVERSE( image[*,*,*,ni],  3,  /Overwrite)
			ENDIF
		ENDFOR
	ENDIF
ENDELSE

error_id = 0
RETURN, REFORM(image)

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
;+
;
;NAME:
;   Read_Raw_v2
;
; -- PROPOSITO
;
;   Lee datos en formato binario (RAW)
;
; -- EJEMPLO DE LLAMADA EN LINEA DE COMANDOS IDL:
;
;
; datos = Read_raw($
;    FILE=file, $
;    OFFSET=offset,
;    DIMENSIONS=dimensions,
;    SWAP=swap,
;    TYPE=type, $
;    REVERSE_XYZ=reverse_xyz,
;    ERROR_ID = error_id,
;    ERROR_STR = error_str)
;
; -- ENTRADAS:
;
;FILE
;   Archivo de lectura de los datos. tipo STRING
;
;OFFSET
;   offset al comienzo del fichero, en número de bytes
;
;DIMENSIONS
;   array de valores con las dimensiones de los datos (hasta cuatro dimensiones, X,Y,Z,T)
;
;SWAP
;   Opción para leer en little endian
;
;MATCH_SIZE
;   Solo lee los datos si encaja con las dimensiones dadas
;
;TYPE
;   variable de tipo STRING o bien INTEGER que especifica el tipo de datos
;   Admite los tipos IDL:
;
;	1 ó 'BYTE'     byte (8 bits)
;	2 ó 'INT':     Entero corto (2 bytes)
;	2 ó 'INTEGER': Entero corto (2 bytes)
;	3 ó 'LONG':    Entero largo (4 bytes)
;	4 ó 'FLOAT':   Punto flotante de 4 bytes
;	4 ó 'FLOATING':Punto flotante de 4 bytes
;	5 ó 'DOUBLE':  Punto flotante de 8 bytes
;	12 ó 'UINT':    Entero corto sin signo (2 bytes)
;	13 ó 'ULONG':   Entero largo sin signo (4 bytes)
;	14 ó 'LONG64':  Entero extralargo (8 bytes)
;	15 ó 'ULONG64': Entero extralargo sin signo (8 bytes)
;
;-


FUNCTION Read_raw_v2, FILE=file, OFFSET=offset, DIMENSIONS=dimensions, SWAP=swap, TYPE=type, $
    REVERSE_ARR=reverse_arr, ERROR_ID = error_id, ERROR_STR = error_str, MATCH_SIZE=match_size


error_id=0
error_str = ''
IF N_ELEMENTS(type) EQ 0 THEN type2 = FIX(4) ELSE type2 = type[0]

opt_match_size = KEYWORD_SET(match_size)

arr_valid_types = [1,2,3,14,4,5,12,13,15L]
arr_bytestype   = [1,2, 4, 8, 4, 8, 2, 4, 8l] ; bytes per data

IF SIZE(type2, /TYPE) EQ 7 THEN BEGIN ; ?type definido como string
	CASE STRUPCASE(STRTRIM(type2,2)) OF
	'BYTE': 	imag_type = 1
	'INT':  	imag_type = 2
	'INTEGER':  imag_type = 2
	'LONG':  	imag_type = 3
	'FLOAT':  	imag_type = 4
	'FLOATING': imag_type = 4
	'DOUBLE':  	imag_type = 5
	'UINT':  	imag_type = 12
	'UINTEGER': imag_type = 12
	'ULONG':  	imag_type = 13
	'LONG64':  	imag_type = 14
	'ULONG64':  imag_type = 15
	ELSE : BEGIN
		error_id = 1 & error_str = 'Not valid type:' + type2
		RETURN, -1
	END
	ENDCASE
ENDIF ELSE BEGIN
	imag_type=FIX(type2[0])
	IF MAX(imag_type EQ arr_valid_types) NE 1 THEN BEGIN
		error_id = 1 & error_str = 'Not valid type:' + STRTRIM(STRING(imag_type),2)
		RETURN, -1
	ENDIF
ENDELSE

IF N_ELEMENTS(reverse_arr) EQ 0 THEN BEGIN
	reverse_arr = BYTARR(N_ELEMENTS(dimensions))
ENDIF ELSE BEGIN
	IF N_ELEMENTS(reverse_arr) NE N_ELEMENTS(dimensions) THEN BEGIN
		error_id = 1 & error_str = 'Not valid inversion array'
		RETURN, -1
	ENDIF
ENDELSE

IF N_ELEMENTS(file) NE 1 THEN RETURN, -1
IF SIZE(file, /TYPE) NE 7 THEN RETURN, -1
IF N_ELEMENTS(dimensions) LE 0 THEN RETURN, -1
IF N_ELEMENTS(dimensions) GT 4 THEN RETURN, -1
IF N_ELEMENTS(dimensions) LT 4 THEN nframes = 1L ELSE nframes = (LONG(dimensions[3]) > 1)
IF N_ELEMENTS(dimensions) LT 3 THEN zsize   = 1L ELSE zsize   = (LONG(dimensions[2]) > 1)
xsize = (LONG(dimensions[0]) > 1)
IF N_ELEMENTS(dimensions) GT 1 THEN $
	ysize = (LONG(dimensions[1]) > 1) ELSE ysize = 1

IF N_ELEMENTS(offset) EQ 0 THEN bytes_offset = 0L ELSE bytes_offset = LONG(offset)

IF FILE_TEST(file[0]) NE 1 THEN BEGIN
	error_str = 'File can not be readed' &   error_id = 2
	RETURN, -1
END

arrayAux  = MAKE_ARRAY(xsize, ysize, TYPE=imag_type)
data = MAKE_ARRAY(xsize, ysize, zsize, nframes, TYPE=imag_type)
data = REFORM(data, xsize, ysize, zsize, nframes)

bytes_per_data = arr_bytestype[(WHERE(imag_type EQ arr_valid_types))[0]]

size_in_bytes  = bytes_per_data*xsize*ysize*zsize*nframes + bytes_offset
; Calculo cual es el tamaño en bytes para asegurarme de que el fichero es del tamaño adecuado

file_size = (FILE_INFO(file[0])).size

IF opt_match_size AND (file_size NE size_in_bytes) THEN BEGIN
	error_id=4 & error_str = 'File size does not match with dimensions' & RETURN, -1
ENDIF
IF size_in_bytes GT file_size THEN BEGIN
	error_id=4 & error_str = 'File size smaller than dimensions' & RETURN, -1
ENDIF
;---------------------------------------------------------------------------------


OPENR, unit, file[0], /GET_LUN, ERROR=err
IF err NE 0 THEN BEGIN
    IF N_ELEMENTS(unit) GT 0 THEN BEGIN & FREE_LUN, unit & CLOSE, unit & ENDIF
    error_id = 1 & error_str = 'Error opening file'
    RETURN, -1
ENDIF
POINT_LUN,  unit,  bytes_offset

FOR ni =0, nframes-1 DO BEGIN
    FOR nz=0, zsize-1 DO BEGIN
       READU, unit, arrayAux
       data[*, *, nz, ni] = arrayAux
    ENDFOR
ENDFOR
FREE_LUN, unit & CLOSE, unit
;----------------------------------------------------------------------------------
IF KEYWORD_SET(swap) THEN data=SWAP_ENDIAN(data)


IF nframes EQ 1 THEN BEGIN
	data = REFORM(data, xsize, ysize, zsize)
	IF  reverse_arr[0] THEN data=REVERSE(data,  1,  /OVERWRITE)
	IF  reverse_arr[1] THEN data=REVERSE(data,  2,  /OVERWRITE)
	IF zsize GT 1 THEN BEGIN
		IF reverse_arr[2] THEN data=REVERSE(data,  3,  /OVERWRITE)
	ENDIF ELSE BEGIN
		data = REFORM(data, xsize, ysize)
	ENDELSE
ENDIF ELSE BEGIN
	IF MAX(reverse_arr) GT 0 THEN BEGIN
		FOR ni=0l, nframes-1 DO BEGIN
			IF  reverse_arr[0] THEN data[*,*,*,ni] = REVERSE(data[*,*,*,ni],  1,  /OVERWRITE)
			IF  reverse_arr[1] THEN data[*,*,*,ni] = REVERSE(data[*,*,*,ni],  2,  /OVERWRITE)
			IF zsize GT 1 THEN BEGIN
				IF reverse_arr[2] THEN  data[*,*,*,ni] = REVERSE( data[*,*,*,ni],  3,  /OVERWRITE)
			ENDIF
		ENDFOR
	ENDIF
ENDELSE
IF reverse_arr[3] THEN BEGIN
	data = REVERSE( data,  4,  /OVERWRITE)
ENDIF

error_id = 0
RETURN, data

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION ReadArray_numbers, array, CHAR_SEP=char_sep, CHAR_BEGIN=char_begin, CHAR_END=char_end,$
INTEGER=integer

IF N_ELEMENTS(array) NE 1 THEN RETURN, -1
IF N_ELEMENTS(char_sep) EQ 0 THEN str_sep = ',' ELSE str_sep = STRING(char_sep[0])
IF N_ELEMENTS(char_begin) EQ 0 THEN str_begin = '[' ELSE str_begin = STRING(char_begin[0])
IF N_ELEMENTS(char_end) EQ 0 THEN str_end = ']' ELSE str_end = STRING(char_sep[0])

;---------------------------------------------------------------------------------
	array_temp = STRTRIM(array,2)
    array_temp = STRSPLIT(array_temp, str_begin, /EXTRACT, /PRESERVE_NULL)
    IF N_ELEMENTS(array_temp) NE 2 THEN RETURN, -1
    array_temp = array_temp[1]
    array_temp = STRSPLIT(array_temp, str_end, /EXTRACT, /PRESERVE_NULL)
    IF N_ELEMENTS(array_temp) NE 2 THEN RETURN, -1
    array_temp = array_temp[0]
    array_temp = STRSPLIT(array_temp, str_sep, /EXTRACT, /PRESERVE_NULL)

    IF KEYWORD_SET(integer) THEN BEGIN
    	result = FIX(array_temp)
    ENDIF ELSE BEGIN
    	result = FLOAT(array_temp)
    ENDELSE
    ;---------------------------------------------------------------------------------
	RETURN, result

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION ReadFile_ASCII, str_file, REMOVE_BLANKS=remove_blanks, ERROR_ID=error_id, $
    MAX_LINES=max_lines, MAX_BYTES=max_bytes

; lee el contenido de un archivo de texto
; limitado a 1000 líneas, o max_lines

error_id=0

IF N_ELEMENTS(str_file) NE 1  THEN BEGIN & error_id = 1 & RETURN, -1 & ENDIF
IF SIZE(str_file, /TYPE) NE 7 THEN BEGIN & error_id = 2 & RETURN, -1 & ENDIF

IF KEYWORD_SET(remove_blanks) THEN opt_removeblanks=1 ELSE opt_removeblanks=0

IF N_ELEMENTS(str_file) NE 1 THEN BEGIN & error_id = 3 & RETURN, -1 & ENDIF
IF N_ELEMENTS(max_lines) NE 1 THEN maximum_lines = 10000 ELSE maximum_lines = LONG(max_lines[0])
IF N_ELEMENTS(max_bytes) NE 1 THEN maximum_bytes = 1000000l ELSE maximum_bytes = LONG(max_bytes[0])
IF FILE_TEST(str_file, /READ, /REGULAR) NE 1 THEN BEGIN & error_id = 1 & RETURN, 1 & ENDIF

OPENR, unit, str_file, /GET_LUN, ERROR= error
IF error NE 0 THEN BEGIN
    IF N_ELEMENTS(unit) GT 0 THEN BEGIN & FREE_LUN, unit & CLOSE, unit & ENDIF
    error_id = 4 & RETURN, -1
ENDIF

st_unit = FSTAT(unit)
IF st_unit.size GE maximum_bytes THEN BEGIN
    FREE_LUN, unit & CLOSE, unit &  error_id = 5 & RETURN, -1
ENDIF

str_line = ''
i=0L
strarr_lines = STRARR(1,maximum_lines)
WHILE (EOF(unit) EQ 0) AND (i LE maximum_lines) DO BEGIN
	READF, unit, str_line
	str_len = STRLEN(str_line)
	IF opt_removeblanks EQ 0 THEN BEGIN
		strarr_lines[i] = str_line
		i++
    ENDIF ELSE BEGIN
		IF str_len GT 0 THEN BEGIN
			strarr_lines[i] = str_line
	        i++
		ENDIF
    ENDELSE
ENDWHILE
FREE_LUN, unit & CLOSE, unit

IF i EQ 0 THEN RETURN, ''
strarr_lines = strarr_lines[0, 0:i-1]

RETURN, strarr_lines
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION ReadStrarr_obtainvalues, strarr_lines, LOWCASE=lowcase, UPCASE=upcase, $
    STR_SEPARATOR = str_separator, ERROR_ID=error_id, ERROR_STR=error_str

; Para obtener una serie de parejas de strings ('identificador', variable') y eliminar
; comentarios de un array de strings (normalmente leido de disco)

error_id = 0
error_str = ''

str_comment = ';'
IF N_ELEMENTS(str_separator) EQ 0 THEN str_separator = ':='

n_lines = N_ELEMENTS(strarr_lines)

FOR i=0, n_lines-1 DO BEGIN ;Erase comments
    str_linesplit  = STRSPLIT(strarr_lines[i], str_comment, /EXTRACT, /PRESERVE_NULL)
    strarr_lines[i] = str_linesplit[0]
ENDFOR

;-------------------------------------------
j=0l
strarr_values = STRARR(2)

FOR i=0l, n_lines-1 DO BEGIN    ;Erase blank lines
    pos_sep = STRPOS(strarr_lines[i], str_separator)
    IF pos_sep NE -1 THEN BEGIN
       str_linesplit = STRSPLIT(strarr_lines[i], str_separator, /EXTRACT, /REGEX)
       CASE  N_ELEMENTS(str_linesplit) OF
         1 :  str_linesplit=[str_linesplit, '']
         2 :
         ELSE : BEGIN &   error_id = 1 &  error_str = 'Not valid config format in parameter: "' + STRTRIM(strarr_lines[i],2)  + '"' &  RETURN, -1 & END
       ENDCASE
       strarr_values[*,j] = str_linesplit
       strarr_values   = [[strarr_values],['','']]
       j=j+1
    ENDIF ELSE BEGIN
    	error_id = 2 & error_str = 'Not valid config format in parameter: "' + STRTRIM(strarr_lines[i],2) + '"' & RETURN, -1
	ENDELSE
ENDFOR
;-------------------------------------------
n_lines = N_ELEMENTS(strarr_values[1,*])-1
IF n_lines GT 0 THEN $
    strarr_values = strarr_values[*, 0:n_lines-1] $
ELSE strarr_values = ''
;-------------------------------------------

strarr_values = STRTRIM(strarr_values,2)    ; Quita espacios en blanco (principio y final)

IF KEYWORD_SET(lowcase) THEN $
    strarr_values = STRLOWCASE(strarr_values)  ; Pasa a minusculas
IF KEYWORD_SET(upcase) THEN $
    strarr_values = STRUPCASE(strarr_values)   ; Pasa a minusculas

RETURN, strarr_values

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>


FUNCTION RelativeContrastEnhancement, data, FRAME_INJECTION=frame_injection, IDX_BEFORE = idx_before, IDX_AFTER=idx_after, $
	MIN_THRESHOLD=min_threshold, MAX_THRESHOLD=max_threshold, ABSOLUTEMINMAX=absoluteminmax

;  a partir de una secuencia de imágenes dinámicas retorna la imagen de contrast enhancement,
; CON RESPECTO A LA PRIMERA IMAGEN O CONJUNTO DE PRIMERAS IMAGENES

; solo contraste positivo...

; ABSOLUTEMINMAX -> calcula el contraste absoluto entre mínimo y máximo

opt_absoluteminmax= KEYWORD_SET(absoluteminmax)

ndim = SIZE(data, /N_DIMENSIONS)

opt_min_threshold = N_ELEMENTS(min_threshold) NE 0
opt_max_threshold = N_ELEMENTS(max_threshold) NE 0


IF N_ELEMENTS(frame_injection) THEN BEGIN
	idx_before = INDGEN(frame_injection+1)
ENDIF ELSE BEGIN
	IF N_ELEMENTS(idx_before) EQ 0 THEN idx_before = 0
ENDELSE

data_temp = data

IF N_ELEMENTS(idx_after) EQ 0 THEN opt_after=0 ELSE opt_after = 1

xsize = (SIZE(data, /DIMENSIONS))[0]
ysize = (SIZE(data, /DIMENSIONS))[1]

CASE ndim OF
3 : BEGIN ; Imagen, un slice
	IF opt_after EQ 1 THEN idx_after  = INDGEN((SIZE(data, /DIMENSIONS))[2])
	im_rce = FLTARR(xsize, ysize)

	CASE opt_absoluteminmax OF
	0 : BEGIN
		FOR i=0l, xsize-1 DO BEGIN
			FOR j=0l, ysize-1 DO BEGIN
				arr_data    = REFORM(data_temp[i,j,*])
				mean_before = MEAN(arr_data[idx_before])
				IF opt_min_threshold THEN mean_before = mean_before > min_threshold[0]
				IF opt_after THEN max_value = MAX(arr_data[idx_after]) ELSE $
					max_value   = MAX(arr_data)
				IF opt_max_threshold THEN max_value = max_value < max_threshold[0]
				im_rce[i,j] = max_value/mean_before
			ENDFOR
		ENDFOR
	END
	1 : BEGIN
		FOR i=0l, xsize-1 DO BEGIN
			FOR j=0l, ysize-1 DO BEGIN
				arr_data    = REFORM(data_temp[i,j,*])
				min_value   = MIN(arr_data, MAX=max_Value)
				IF opt_min_threshold THEN min_value  = min_value > min_threshold[0]
				IF opt_max_threshold THEN max_value = max_value < max_threshold[0]
				im_rce[i,j] = max_value/min_value
			ENDFOR
		ENDFOR
	END

	ENDCASE
END
4 : BEGIN ; Imagen 3D ( varios cortes)
	n_slices = (SIZE(data, /DIMENSIONS))[2]
	IF opt_after EQ 1 THEN idx_after  = INDGEN((SIZE(data, /DIMENSIONS))[3])
	im_rce = FLTARR(xsize, ysize, n_slices)

	CASE opt_absoluteminmax OF
	0 : BEGIN
		FOR i=0l, xsize-1 DO BEGIN
			FOR j=0l, ysize-1 DO BEGIN
				FOR k=0l, n_slices-1 DO BEGIN
					arr_data = REFORM(data_temp[i,j,k,*])
					mean_before = MEAN(arr_data[idx_before])
					IF opt_after THEN max_value = MAX(arr_data[idx_after]) ELSE $
						max_value   = MAX(arr_data)
					im_rce[i,j,k] = max_value/mean_before
				ENDFOR
			ENDFOR
		ENDFOR
	END
	1 : BEGIN
		FOR i=0l, xsize-1 DO BEGIN
			FOR j=0l, ysize-1 DO BEGIN
				FOR k=0l, n_slices-1 DO BEGIN
					arr_data = REFORM(data_temp[i,j,k,*])
					min_value   = MIN(arr_data, MAX=max_Value)
					im_rce[i,j,k] = max_value/min_value
				ENDFOR
			ENDFOR
		ENDFOR
	END
	ENDCASE

END
ELSE : RETURN, -1
ENDCASE

RETURN ,im_rce

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Remove_Comments, array, STR_COMMENTS=str_comments

n_lines = N_ELEMENTS(array)

array_sal = STRARR(1, n_lines)


n=0l
FOR i=0l, n_lines-1 DO BEGIN
	pos = (STRPOS(array[i], str_comments))[0]
	IF pos NE -1 THEN BEGIN
		array_temp = STRMID(array[i],0,pos)
	ENDIF ELSE BEGIN
		array_temp = array[i]
	ENDELSE
	IF STRLEN(array_temp) GT 0 THEN BEGIN
		array_sal[n]=array_temp
		n++
	ENDIF
ENDFOR
array_sal = array_sal[0,0:n-1]

RETURN, array_sal


END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Restore_GraphicObjectsData, FILE=file

	IF FILE_TEST(file) NE 1 THEN RETURN, -1

	RESTORE, file
	IF N_ELEMENTS(stobjects) EQ 0 THEN RETURN,-1

	;stobjects = {configplots  :configplots,   $
	;		configaxis    :configaxis,    $
	;		configaxistext:configaxistext,$
	;		configfont    :configfont,    $
	;		configsymbol  :configsymbol,  $
	;		configtext    :configtext}

	RETURN, stobjects

	;stplot1 =  PlotStruct_GetValues(stobjects.configplots[0])
	;stplot2 =  PlotStruct_GetValues(stobjects.configplots[1])

	;RETURN, [stplot1, stplot2]

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION ROImap_GetPoint, ST_ROI=st_roi, point, POINT_IMF=point_imf


val_idx = (*st_roi.im_idxs)[point[0],point[1]]

IF val_idx NE -1 THEN BEGIN

	; Remember, This function works in IDL standard coordinates (origin of image (x=0,y=0) is in the lower-left corner)

	point_diff_x_lr = 1.0*(point[0]-st_roi.pr_ini[0])/st_roi.np_lr
	point_diff_y_lr = 1.0*(point[1]-st_roi.pr_ini[1])/st_roi.np_lr

	IF (point_diff_x_lr GE 0) AND (point_diff_y_lr GE 0) AND $
		(point_diff_x_lr LE (st_roi.roi_sz[0]/st_roi.np_lr)) AND $
		(point_diff_y_lr LE (st_roi.roi_sz[1]/st_roi.np_lr)) THEN BEGIN

		ppx = 0 > FIX(point_diff_x_lr) < ((st_roi.roi_sz[0]/st_roi.np_lr)-1)
		ppy = 0 > FIX(point_diff_y_lr) < ((st_roi.roi_sz[1]/st_roi.np_lr)-1)

		point_imf = [ppx, (st_roi.roi_sz[1]/st_roi.np_lr)-1-ppy]

	ENDIF
ENDIF

RETURN, val_idx

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION SetWidget_maxsizes, arr_ids, XSIZE=xsize, YSIZE=ysize

opt_xsize = KEYWORD_SET(xsize)
opt_ysize = KEYWORD_SET(ysize)

n_widgets = N_ELEMENTS(arr_ids)
arr_sizes = LONARR(n_widgets)

FOR I=0l, n_widgets-1 DO BEGIN
	IF opt_xsize THEN BEGIN
		arr_sizes[i] = (WIDGET_INFO(arr_ids[i], /GEOMETRY)).scr_xsize
	ENDIF ELSE BEGIN
		arr_sizes[i] = (WIDGET_INFO(arr_ids[i], /GEOMETRY)).scr_ysize
	ENDELSE
ENDFOR
max_size = MAX(arr_sizes)

FOR I=0l, n_widgets-1 DO BEGIN
	IF opt_xsize THEN BEGIN
		WIDGET_CONTROL, arr_ids[i], SCR_XSIZE=max_size
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, arr_ids[i], SCR_YSIZE=max_size
	ENDELSE
ENDFOR

RETURN, max_size

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION SetWidget_RelativePosition, base, ID_PARENT=id_parent

	arr_positions = [10,14,9,13,2,8,0,6,1,7,4]
	np=0L
	WHILE  N_ELEMENTS(position) NE 2 DO BEGIN
		position = GetWidget_RelativePosition(base, ID_PARENT=id_parent, POSITION=arr_positions[np], /NOOUTSIDE)
		np++
		IF np GT 10 THEN position = [0,0l]
	ENDWHILE
	WIDGET_CONTROL, base, XOFFSET=position[0], YOFFSET=position[1]

	;  9|0   1   2|10
	;   |         |
	; 11|3   4   5|12
	;   |         |
	; 13|6   7   8|14
	;    ---------

RETURN, 1

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>


FUNCTION Sigmoid_curve, arr_t, DELAY_T=delay_t, ALPHA=alpha, TAU=tau

;---------------------------------------------------------------------------------
;M. C. Schabel, and D. L. Parker,
;Uncertainty and bias in contrast concentration measurements using spoiled gradient echo pulse sequences,
;Physics in Medicine and Biology, vol. 53, no. 9, pp. 2345-2373, May 7, 2008.
;---------------------------------------------------------------------------------
; (EQUATION A.3)

; Double precission
alphad = DOUBLE(alpha[0])
taud   = DOUBLE(tau[0])
delayd = DOUBLE(delay_t[0])


pos_1 = WHERE(arr_t GT 0, ct1)
pos_2 = WHERE(arr_t LE 0, ct2)

arr_result = DBLARR(N_ELEMENTS(arr_t))


IF ct1 GT 0 THEN BEGIN
	arr_t_1    = arr_t[pos_1]

	IF taud EQ 0 THEN BEGIN
		igamma_result = 1d
	ENDIF ELSE BEGIN
		par_aux = (1.0d/taud - 1.0d/delayd)*arr_t_1
		IF MIN(par_aux) LT 0 THEN BEGIN
			RETURN, -1
		ENDIF
		igamma_result =IGAMMA(alphad, par_aux)
	ENDELSE


	;igamma_result = IGAMMA(alphad, par_aux)
	;IF N_ELEMENTS(igamma_result) NE N_ELEMENTS(par_aux) THEN BEGIN
	;	PRINT, 'Error'
	;	RETURN, -1
	;ENDIF

	arr_result_1 = (delayd/((delayd-taud)*GAMMA(alphad))) * $
			(arr_t_1/delayd)^(-taud/(delayd-taud)) * $
			(EXP(-arr_t_1/delayd)) * $
			igamma_result*GAMMA(alphad)


	arr_result[pos_1] = arr_result_1
ENDIF
IF ct2 GT 0 THEN BEGIN
	arr_result[pos_2] = 0
ENDIF

RETURN, arr_result



END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION StandardError_divide, MEANS=means, SD=sd

; standard deviation of multiplication of two variables with known STANDARD DEVIATION AND MEAN

std_error_total = SQRT((sd[0]/means[0])^2 + (sd[1]/means[1])^2)*(means[0]/means[1])


RETURN, std_error_total

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION StandardError_multiply, MEANS=means, SD=sd

; standard deviation of multiplication of two variables with known STANDARD DEVIATION AND MEAN

std_error_total = SQRT((sd[0]/means[0])^2 + (sd[1]/means[1])^2)*(means[0]*means[1])


RETURN, std_error_total

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION Translate_data2strings, data, N_DECIMALS=n_decimals, REVERSE_LINES=reverse_lines, INTEGER=integer,$
	MIN_VALUE=min_value, MAX_VALUE=max_value, STR_NOVALID=str_novalid,$
	HEADER=header, COL_TITLES=col_titles, MIN_COL_LENGTH=min_col_length


IF N_ELEMENTS(data) eq 0 THEN	RETURN,-1

opt_minval = N_ELEMENTS(min_value) NE 0
opt_maxval = N_ELEMENTS(max_value) NE 0

CASE SIZE(data, /N_DIMENSIONS) OF
2 : BEGIN
	n_col = (SIZE(data, /DIMENSIONS))[0]
	n_fil = (SIZE(data, /DIMENSIONS))[1]
END
1 : BEGIN
	n_col = 1l
	n_fil = N_ELEMENTS(data)
END
ELSE : RETURN, -1
ENDCASE

IF opt_minval EQ 1 	AND n_col GT 1 THEN BEGIN
	IF N_ELEMENTS(min_Value) EQ 1 THEN min_value=REPLICATE(min_value, n_col)
	IF N_ELEMENTS(min_value) NE n_col THEN RETURN,-1
ENDIF
IF opt_maxval EQ 1 	AND n_col GT 1 THEN BEGIN
	IF N_ELEMENTS(max_Value) EQ 1 THEN max_value=REPLICATE(max_value, n_col)
	IF N_ELEMENTS(max_value) NE n_col THEN RETURN,-1
ENDIF

arr_strings = STRARR(n_col, n_fil)

IF NOT KEYWORD_SET(integer) THEN $
	st_format = '(F30.' + STRTRIM(STRING(n_decimals),2) + ")"
; los datos seran enteros si no,


FOR j=0l, n_fil-1 DO BEGIN
	FOR i=0l, n_col-1 DO BEGIN
		val = data[i+j*n_col]
		str_number = STRTRIM(STRING(val, FORMAT=st_format),2)
		IF opt_minval THEN BEGIN
			IF val LT min_value[i] THEN str_number = str_novalid
		ENDIF
		IF opt_maxval THEN BEGIN
			IF val GT max_value[i] THEN str_number = str_novalid
		ENDIF
		arr_strings[i+j*n_col] = str_number
	ENDFOR
ENDFOR

arr_max_length = INTARR(n_col)
FOR i=0l, n_col-1 DO BEGIN
	arr_max_length[i] = MAX(STRLEN(arr_strings[i,*]))
ENDFOR
arr_max_length[*] = MAX(arr_max_length)

IF N_ELEMENTS(min_col_length) NE 0 THEN arr_max_length[*] = arr_max_length[*] > min_col_length


;max_length = MAX(STRLEN(arr_strings))

FOR J=0L, n_fil-1 DO BEGIN
	FOR i=0L, n_col-1 DO BEGIN
		IF STRLEN(arr_strings[j*n_col+i]) LT arr_max_length[i] THEN BEGIN
			n_spaces = arr_max_length[i] - STRLEN(arr_strings[j*n_col+i])
			IF n_spaces EQ 0 THEN PRINT, 'EEEEE'
			arr_strings[j*n_col+i] = STRING(BYTARR(n_spaces)+32b) + arr_strings[j*n_col+i]
		ENDIF
		IF i NE n_col-1 THEN arr_strings[i,j]+=' ' ;STRING(9B)
	ENDFOR
ENDFOR

IF n_col GT 1 THEN $
	arr_strings = STRJOIN(arr_strings)

IF KEYWORD_SET(reverse_lines) THEN BEGIN
	n_dim = SIZE(arr_strings, /N_DIMENSIONS)
	arr_strings = REVERSE(arr_strings, n_dim)
ENDIF

IF N_ELEMENTS(col_titles) NE 0 THEN BEGIN
	n_col_tit = N_ELEMENTS(col_titles)
	strarr_coltitles = STRARR(n_col_tit)
	FOR i=0L, n_col_tit-1 DO BEGIN
		IF STRLEN(col_titles[i]) LT arr_max_length[i] THEN BEGIN
			n_spaces = arr_max_length[i] - STRLEN(col_titles[i])
			strarr_coltitles[i] = STRING(BYTARR(n_spaces)+32b) + col_titles[i]
		ENDIF ELSE BEGIN
			strarr_coltitles[i] = STRMID(col_titles[i],0,arr_max_length[i])
		ENDELSE
		IF i NE n_col-1 THEN strarr_coltitles[i]+=' ' ;STRING(9B)
	ENDFOR
	str_coltitles = STRJOIN(strarr_coltitles)
	arr_strings = [str_coltitles, arr_strings]
ENDIF


IF N_ELEMENTS(header) GT 0 THEN BEGIN
	arr_strings = [header, arr_strings]
ENDIF

RETURN, arr_strings

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Translate_number2string, NUMBER=number, LENGTH=length

; Devuelve un numero entero pasado a string, precedido de ceros hasta llegar a length

num = FIX(number)
str_num = STRTRIM(STRING(num),2)
len_num = STRLEN(str_num)

str_zeros = ''
diff = (length-len_num) > 0
FOR i=0, diff-1 DO BEGIN
	str_zeros = str_zeros + '0'
ENDFOR
str_num = str_zeros + str_num
RETURN, str_num

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Translate_numberarray2string, ARRAY=array, FLOAT_N=float_n, SEP=sep, BRACKETS=brackets

; Devuelve un array de un solo elemento, de una serie de numeros
;float_n : Numero de elementos de punto flotante despues del punto

n_elem = N_ELEMENTS(array)
IF N_ELEMENTS(sep) EQ 0 THEN sep = ' '
IF N_ELEMENTS(float_n) EQ 0 THEN float_length=5 ELSE float_length=FIX(float_n[0])

IF n_elem LE 0 THEN RETURN, ''
type_number = SIZE(array, /TYPE)
IF (type_number EQ 4 OR type_number EQ 5) THEN BEGIN
	IF float_length GT 0 THEN BEGIN
		ff = '(F' + STRTRIM(STRING(float_length+20),2) + '.' + STRTRIM(STRING(float_length),2) + ')'
		number_str = STRTRIM(STRING(array, FORMAT=ff),2)
	ENDIF ELSE BEGIN
		number_str = STRTRIM(STRING(FIX(array)),2)
	ENDELSE
ENDIF ELSE BEGIN
	number_str = STRTRIM(STRING(array),2)
ENDELSE

number_ss = ''
FOR i=0, n_elem-1 DO BEGIN
	number_ss = number_ss + number_str[i]
	IF i NE n_elem-1 THEN number_ss+=sep
ENDFOR

IF KEYWORD_SET(brackets) THEN number_ss = '['+number_ss + ']'


RETURN, number_ss
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
PRO Undefine, var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, $
	var11, var12, var13, var14, var15, var16, var17, var18, var19, var20

; Procedimiento de borrado de variables liberación de memoria
; Limitada a 20 variables como máximo. No libera la memoria de variables de tipo objeto y puntero

;n_t = N_PARAMS()

IF N_ELEMENTS(var1) NE 0 THEN temp = TEMPORARY(var1) ;ELSE RETURN
IF N_ELEMENTS(var2) NE 0 THEN temp = TEMPORARY(var2) ;ELSE RETURN
IF N_ELEMENTS(var3) NE 0 THEN temp = TEMPORARY(var3) ;ELSE RETURN
IF N_ELEMENTS(var4) NE 0 THEN temp = TEMPORARY(var4) ;ELSE RETURN
IF N_ELEMENTS(var5) NE 0 THEN temp = TEMPORARY(var5) ;ELSE RETURN
IF N_ELEMENTS(var6) NE 0 THEN temp = TEMPORARY(var6) ;ELSE RETURN
IF N_ELEMENTS(var7) NE 0 THEN temp = TEMPORARY(var7) ;ELSE RETURN
IF N_ELEMENTS(var8) NE 0 THEN temp = TEMPORARY(var8) ;ELSE RETURN

IF N_ELEMENTS(var9)  NE 0 THEN temp = TEMPORARY(var9)  ;ELSE RETURN
IF N_ELEMENTS(var10) NE 0 THEN temp = TEMPORARY(var10) ;ELSE RETURN
IF N_ELEMENTS(var11) NE 0 THEN temp = TEMPORARY(var11) ;ELSE RETURN
IF N_ELEMENTS(var12) NE 0 THEN temp = TEMPORARY(var12) ;ELSE RETURN
IF N_ELEMENTS(var13) NE 0 THEN temp = TEMPORARY(var13) ;ELSE RETURN
IF N_ELEMENTS(var14) NE 0 THEN temp = TEMPORARY(var14) ;ELSE RETURN
IF N_ELEMENTS(var15) NE 0 THEN temp = TEMPORARY(var15) ;ELSE RETURN
IF N_ELEMENTS(var16) NE 0 THEN temp = TEMPORARY(var16) ;ELSE RETURN
IF N_ELEMENTS(var17) NE 0 THEN temp = TEMPORARY(var17) ;ELSE RETURN
IF N_ELEMENTS(var18) NE 0 THEN temp = TEMPORARY(var18) ;ELSE RETURN
IF N_ELEMENTS(var19) NE 0 THEN temp = TEMPORARY(var19) ;ELSE RETURN
IF N_ELEMENTS(var20) NE 0 THEN temp = TEMPORARY(var20) ;ELSE RETURN

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

FUNCTION UpdateROIdata, IDX_POINTS=idx_points, ST_ROI=st_roi, IMAGESIZE=imagesize, ERROR_ID=error_id, ERROR_STR=error_str, SLICEZ=slicez


error_id  = 0
error_str = ''

;-----------------------------------------------------------------------------------------------
; valido para todo tipo de ROIS a partir de st.results.idx_points
idx_points_lr = IdxPointsLr_From_Indexes(IDX_POINTS=idx_points, $
	ROI_RESOLUTION=st_roi.np_lr, MIN_P=min_p, MAX_P=max_p, DIMENSIONS=imagesize, $
	MASK_LR=mask_lr, MASK_HR=mask_hr, ERROR_ID=error_id)
;-----------------------------------------------------------------------------------------------
IF error_id EQ 2 THEN BEGIN
	error_str = 'No points in ROI, please repeat' & RETURN,-1 & ENDIF
IF error_id NE 0 THEN BEGIN
	error_str = 'Error selection ROI' & RETURN,-1 & ENDIF
;-----------------------------------------------------------------------------------------------

st_roi.pr_ini         = min_p
st_roi.roi_sz         = max_p-min_p+1
st_roi.slicez         = slicez

*st_roi.mask_lr       = mask_lr
*st_roi.idx_points_lr = idx_points_lr
*st_roi.idx_points    = idx_points

;-----------------------------------------------------------------------------------------------
*st_roi.im_idxs = ImIdxs_From_Indexes(MASK_LR=mask_lr,DIMENSIONS=imagesize,$
		ROI_RESOLUTION=st_roi.np_lr, MIN_P = st_roi.pr_ini)
;-----------------------------------------------------------------------------------------------

RETURN, st_roi

END

;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>

PRO view, imagen, wind, CONTRAST=contrast, NOSCALE=noscale

;; Visualiza correctamente una imagen en pantalla
;; En la ventana indicada (cero por defecto)

IF N_ELEMENTS(wind) EQ 0  THEN wind=0	; Ventana por defecto
IF (SIZE(imagen))[0] NE 2 THEN BEGIN
	PRINT, 'datos incorrectos'
	RETURN
ENDIF

IF NOT KEYWORD_SET(noscale) THEN noscal=0 ELSE noscal=1

IF N_ELEMENTS(contrast) NE 2 THEN ct=[0,1d] ELSE BEGIN
	ct = [DOUBLE(contrast[0]), DOUBLE(contrast[1])]
ENDELSE

tam = SIZE(imagen)

lg_win = tam[1]
wd_win = tam[2]

WINDOW, wind,  XSIZE=lg_win, YSIZE=wd_win
IF noscal EQ 0 THEN $
	TVSCL, BYTSCL(imagen, MIN=MIN(imagen)*ct[0], MAX=MAX(imagen)*ct[1]), ORDER=!order $
ELSE $
	TV, imagen, ORDER=!order

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION view_mosaic, arr_images, N_COL=n_col, DIV=div, WINDOW_NUMBER=window_number, NONEW=nonew, $
	YPOS=ypos,XPOS=xpos, BACKGROUND=background, DIFF_IM=diff_im, BYTSCL_EQUAL=bytscl_equal,$
	PALETTE=palette


; MIRA UN MOSAICO DE IMAGENES DEL MISMO TAMAÑO

IF N_ELEMENTS(window_number) EQ 0 THEN window_number = 0
IF KEYWORD_SET(nonew) THEN opt_nonew = 1 ELSE opt_nonew = 0

;-------------------------------------
size_im = SIZE(arr_images, /DIMENSIONS)
xsize = size_im[0]
ysize = size_im[1]
n_images = size_im[2] ; blanco y negro
;-------------------------------------

IF N_ELEMENTS(xpos)  EQ 0 THEN xpos=20
IF N_ELEMENTS(ypos)  EQ 0 THEN ypos=20
IF N_ELEMENTS(div)   EQ 0 THEN divv = 4.0  ELSE divv = FLOAT(div)
IF N_ELEMENTS(n_col) EQ 0 THEN n_col = CEIL(SQRT(n_images))

n_row = CEIL(n_images*1.0/n_col)
;-------------------------------------

IF N_ELEMENTS(background) EQ 0 THEN BEGIN
	!P.BACKGROUND = '804020'x
ENDIF ELSE !P.BACKGROUND = background

tam_x  = xsize/divv & tam_y  = ysize/divv

IF N_ELEMENTS(diff_im) EQ 0 THEN diff_im = tam_X/8.0

IF opt_nonew eq 0 THEN BEGIN
	WINDOW, window_number, XSIZE=tam_x*n_col+ diff_im*(n_col+1), YSIZE=tam_y*n_row + diff_im*(n_row+1), $
		XPOS=xpos, YPOS=ypos
	DEVICE, DECOMPOSED=1
	IF N_ELEMENTS(background) EQ 0 THEN BEGIN
	!P.BACKGROUND = '804020'x
	ENDIF ELSE !P.BACKGROUND = background
	ERASE
	DEVICE, DECOMPOSED=0
ENDIF ELSE BEGIN
	WSET, window_number
ENDELSE

IF KEYWORD_SET(bytscl_equal) THEN BEGIN
	max_bytscl = MAX(arr_images)
ENDIF

IF N_ELEMENTS(palette) GT 0 THEN BEGIN
	TVLCT,rb, gb, bg, /GET
	TVLCT, palette[*,0],palette[*,1],palette[*,2]
ENDIF

FOR i=0, n_images-1 DO BEGIN
	im = arr_images[*,*,i]
	idx_x = i MOD n_col
	idx_y = i/n_col
	im_view = BYTSCL(CONGRID(im, tam_x, tam_y, /INTERP), MAX=max_bytscl)
	TVSCL, im_view, tam_x*idx_x + diff_im*(idx_x+1), tam_y*(n_row-idx_y-1) + diff_im*(n_row-idx_y)
ENDFOR

!P.BACKGROUND = '0'x

IF N_ELEMENTS(palette) GT 0 THEN BEGIN
	TVLCT,rb, gb, bg
ENDIF


RETURN,im_view

WSHOW_ALL

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
PRO viewg, image, ventana, CONTRAST=contrast, INTERP=interp, XSIZE=xsize, YSIZE=ysize, NOSCALE=noscale

; Visualiza una imagen en pantalla a 512*x píxeles

IF N_ELEMENTS(image) EQ 0 THEN RETURN
tam = SIZE(image)
IF tam[0] LT 2 THEN RETURN
true = 0
IF tam[0] EQ 3 THEN BEGIN
	IF MIN(tam[1:3]) GT 2 THEN true = 1 ELSE true = 0
ENDIF
IF true EQ 1 THEN BEGIN
	IF tam[1] EQ 3 THEN imag = change_true(image) ELSE imag = image
ENDIF ELSE imag = REFORM(image)

IF N_ELEMENTS(xsize) EQ 0 THEN tam_v = 512e ELSE tam_v = FIX(xsize[0])

xsize =tam[1]
ysize =tam[2]

tmax = MAX([xsize, ysize])
aspect = FLOAT(ysize)/FLOAT(xsize)

IF N_ELEMENTS(ventana) EQ 0 THEN ventana=0

IF tmax EQ xsize THEN BEGIN
	IF (true EQ 0) THEN BEGIN
	is = congrid(imag, tam_v, FIX(ROUND(tam_v*aspect)), INTERP=interp, /CENTER)
	view,  is, ventana,CONTRAST=contrast, NOSCALE=noscale
ENDIF
IF (true EQ 1) THEN BEGIN
	is = congrid(imag, tam_v,FIX(ROUND(tam_v*aspect)),3,INTERP=interp, /CENTER)
	viewt, is, ventana, NOSCALE=noscale
ENDIF
ENDIF ELSE BEGIN
	IF (true EQ 0) THEN BEGIN
		is = congrid(imag, FIX(ROUND(tam_v/aspect)), tam_v, INTERP=interp, /CENTER)
		view,  is, ventana, CONTRAST=contrast, NOSCALE=noscale
	ENDIF
	IF (true EQ 1) THEN BEGIN
		is = congrid(imag, FIX(ROUND(tam_v/aspect)), FIX(tam_v), 3, INTERP=interp, /CENTER)
		viewt, is, ventana, NOSCALE=noscale
	ENDIF
ENDELSE

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
PRO viewt, imag, wind, NOSCALE=noscale, _EXTRA=extra

;; Visualiza correctamente imágenes de greyscale y color
;;

IF N_ELEMENTS(imag) EQ 0 THEN BEGIN & PRINT, 'Datos no válidos' & RETURN & END
IF NOT KEYWORD_SET(noscale) THEN noscal=0 ELSE noscal=1


tam=SIZE(imag)

IF N_ELEMENTS(wind) EQ 1 THEN BEGIN
	WINDOW, wind   ; Ventana por defecto
ENDIF
IF N_ELEMENTS(wind) GT 1 THEN BEGIN
	PRINT, 'Datos no válidos' & RETURN
ENDIF
IF N_ELEMENTS(wind) EQ 0 THEN BEGIN
	WINDOW, 0  ;
ENDIF

IF tam[0] EQ 2  THEN BEGIN 	; gris
	WINDOW, !D.window,  XSIZE=tam(1), YSIZE=tam(2)
	IF noscal EQ 0 THEN TVSCL, imag, _EXTRA=extra ELSE TV, imag, _EXTRA=extra
	RETURN
ENDIF


IF tam[0] EQ 3  THEN BEGIN     ; color

	IF tam[1] EQ 3  THEN BEGIN  ; TRUE=1
		WINDOW, !D.window,  XSIZE=tam(2), YSIZE=tam(3)
		IF noscal EQ 0 THEN TVSCL, imag, TRUE=1, _EXTRA=extra ELSE TV, imag, TRUE=1, _EXTRA=extra
		RETURN
	ENDIF
	IF tam[2] EQ 3  THEN BEGIN  ; TRUE=2
		WINDOW, !D.window,  XSIZE=tam(1), YSIZE=tam(3)
		IF noscal EQ 0 THEN TVSCL, imag, TRUE=2, _EXTRA=extra ELSE TV, imag, TRUE=3, _EXTRA=extra
		RETURN
	ENDIF
	IF tam[3] EQ 3  THEN BEGIN  ; TRUE=3
		WINDOW, !D.window,  XSIZE=tam(1), YSIZE=tam(2)
		IF noscal EQ 0 THEN TVSCL, imag, TRUE=3, _EXTRA=extra ELSE TV, imag, TRUE=3, _EXTRA=extra
		RETURN
	ENDIF

ENDIF

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
PRO Wdel_all
; Borra todas las ventanas gráficas

j =  0
WHILE j GE 0 DO BEGIN
	i = !D.WINDOW
	j = i
	IF j GE 0 THEN WDELETE, i
ENDWHILE
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION Write_Raw, FILENAME=filename, DATA=data, OVERWRITE=overwrite

opt_overwrite= KEYWORD_SET(overwrite)
IF N_ELEMENTS(data) EQ 0 THEN RETURN, -1
IF N_ELEMENTS(filename) EQ 0 THEN RETURN, -1

IF FILE_TEST(filename) AND opt_overwrite THEN RETURN, -1

IF FILE_TEST(get_path(filename), /DIRECTORY) NE 1 THEN RETURN, -1

OPENW, unit, filename, /GET_LUN, ERROR=err
IF err NE 0 THEN BEGIN & FREE_LUN, unit & RETURN, -1 & ENDIF
WRITEU, unit, data
FREE_LUN, unit
ok = FILE_TEST(filename)
IF ok NE 1 THEN RETURN, -1

RETURN, 1

END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
FUNCTION WriteFile_ASCII, str_file, ARRAY=array, ADD=add, NO_OVERWRITE=no_overwrite, ERROR_ID=error_id

; Escribe un array en un archivo de texto

error_id = 1
IF KEYWORD_SET(add) THEN opt_add=1 ELSE opt_add=0
IF N_ELEMENTS(array) EQ 0 THEN RETURN,-1
type_data = SIZE(array, /TYPE)
IF type_data NE 7 THEN RETURN, -1
num_data  = N_ELEMENTS(array)

IF opt_add EQ 0 THEN BEGIN
    IF KEYWORD_SET(no_overwrite) THEN BEGIN
    	ok = FILE_TEST(str_file, /WRITE, /REGULAR)
    	IF ok EQ 1 THEN RETURN, -1
    ENDIF
    OPENW, unit, str_file, /GET_LUN, ERROR = error
ENDIF ELSE BEGIN
    ok = FILE_TEST(str_file, /WRITE, /REGULAR)
    IF ok NE 1 THEN RETURN, -1
    OPENU, unit, str_file, /GET_LUN, ERROR = error, /APPEND
ENDELSE
IF error NE 0 THEN BEGIN
    RETURN, -1
ENDIF
FOR i=0l, num_data-1 DO BEGIN
    PRINTF, unit, array[i]
ENDFOR
FREE_LUN, unit
CLOSE, unit

error_id = 0
RETURN, 1
END
;<->
;**************************************************************************************************
;**************************************************************************************************
;<+>
PRO Wshow_all
; Muestra todas las ventanas gráficas

DEVICE, WINDOW_STATE=arr_ws
w_valid = WHERE(arr_ws NE 0, ct)
IF ct GT 0 THEN BEGIN
	FOR i=0, ct-1 DO BEGIN
		WSHOW, w_valid[i]
	ENDFOR
ENDIF
END
;<->
;**************************************************************************************************
;**************************************************************************************************

