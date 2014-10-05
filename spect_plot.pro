PRO spect_plot, nt, dt
    dee = 0.05
    nbins = 1000
    ene = FINDGEN(1000)*dee
    ; Open the HDF5 file
    file_id = H5F_OPEN('spectra.h5')
    FOR i = 0, nt, dt DO BEGIN
        dsetname = STRING(i, FORMAT='(I5.5)') + '/spectra'
        dataset_id = H5D_OPEN(file_id, dsetname)
        data = h5D_READ(dataset_id)
        dataspace_id = H5D_GET_SPACE(dataset_id)
        H5S_CLOSE, dataspace_id
        H5D_CLOSE, dataset_id
        p1 = PLOT(ene, data(0,*),$
            /XLOG, /YLOG, $
            FONT_SIZE=16,/OVERPLOT,$
            XTITLE='$E/m_e c^2$',$
            YTITLE='Number of electrons F(E)')
    ENDFOR
    p1.XRANGE=[0.05, 5]
    p1.YRANGE=[1, 1E7]
    H5F_CLOSE, file_id
    p1.save, 'espect.jpg'
;    p1 = PLOT(ts, data(i-1,*),$
;        COLOR=colorlist[i-1],$
;        FONT_SIZE=16,/OVERPLOT,$
;        XRANGE=[xmin, xmax], YRANGE=[ymin,ymax],$
;        XTITLE='$\Omega_{ci}t$', YTITLE="Total Energy",$
;        NAME=namelist[i-1])
;    fname = 'spect.jpg'
;    p1.save, fname
;
END
