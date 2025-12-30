
cd /Users/mikako/work/teaching/2025_ObsTech/LCO/data/NGC3201/NGC3201_B_2025

funpack *.fits.fz

;https://starlink.eao.hawaii.edu/docs/sun55.htx/sun55se2.html#x3-50001
convert

ls -1 *fits


fits2ndf coj0m416-sq36-20251227-0173-e00.fits Bimg1.sdf
fits2ndf coj0m416-sq36-20251227-0173-e91.fits Bimg2.sdf
fits2ndf coj0m416-sq36-20251227-0174-e00.fits Bimg2.sdf
fits2ndf coj0m416-sq36-20251227-0174-e91.fits Bimg3.sdf
fits2ndf coj0m416-sq36-20251227-0175-e00.fits Bimg4.sdf
fits2ndf coj0m416-sq36-20251227-0175-e91.fits Bimg5.sdf
fits2ndf coj0m416-sq36-20251227-0176-e00.fits Bimg6.sdf
fits2ndf coj0m416-sq36-20251227-0176-e91.fits Bimg7.sdf
fits2ndf coj0m416-sq36-20251227-0177-e00.fits Bimg8.sdf
fits2ndf coj0m416-sq36-20251227-0177-e91.fits Bimg9.sdf
fits2ndf coj0m416-sq36-20251227-0178-e00.fits Bimg10.sdf
fits2ndf coj0m416-sq36-20251227-0178-e91.fits Bimg11.sdf


kappa

wcsmosaic
  Bimg*.sdf
  Bimg1.sdf
  Bstack.sdf

  ;This needs jitter correction

/Users/mikako/work/teaching/2024_ObservingTechnique/LCO_2024/LCO/M45_archive_BV/V_long

 wcsmosaic
  M45*.sdf
  M45_Vlong_A.sdf
  [-3064,-6061]
  [3094,2042]
  test.sdf

  ;somehow this works
