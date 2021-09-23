;physical constants
G=6.67430D-8 ; dyne cm^2 g^-2
amu=1.66054d-24 ; g
mperH2=2.8*amu
Msun=1.98847d33 ; g
au=1.495978707d13 ; cm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;read the density cube
nx=long(66) & ny=long(66) & nz=long(66)
dube=READ_BINARY('hmm1_dens.rcube',DATA_TYPE=5,DATA_DIMS=[nx,ny,nz])
;regrid
cdube=CONGRID(dube,nx/2,ny/2,nz/2)
nx=nx/2 & ny=ny/2 & nz=nz/2
dube=cdube

;number of cells
ncells=nx*ny*nz
;physical size of a cell in cm
distinpc=140.0 ; distance in pc
mapsize=130.0 ; map size arcsec
cellsize=mapsize*distinpc*au/float(nz)

;x,y,z coordinates of all cells
z=lindgen(ncells)/(nx*ny)
y=lindgen(ncells)/nx-z*ny
x=lindgen(ncells)-(y+z*ny)*nx

;Mass
mass=total(dube)*cellsize^3*mperH2/Msun
print,format='(A,E12.3,A)', 'Total mass: ', mass, ' Msun'

;gravitational potential cube
Vgrav=dblarr(nx,ny,nz)
;gravitational acceleration 
agrav=dblarr(3,nx,ny,nz)
;directions from one cell to all others 
acomponents=dblarr(3,ncells)

for i=0,ncells-1 do begin
    ;x,y,z coordinates of cell i:
    zc=i/(nx*ny)
    yc=i/nx-zc*ny
    xc=i-(yc+zc*ny)*nx
    ;calculate 1/r_ij
    invdist=1.0/sqrt((float(x-xc))^2+(float(y-yc))^2+(float(z-zc))^2)
    invdist[i]=0.0
    ;V_i = G Sum (m_j/r_ij)
    Vgrav[xc,yc,zc]=total(dube*invdist)
    ;r_ij vectors, this is double array of (3xncells):
    directions=[(x-xc),(y-yc),(z-zc)]
    ;reformat to a double array of (3,ncells)
    ;these are accelerations caused by all other cells: 
    for j=0,2 do begin
      acomponents[j,*]=dube*invdist^3*directions[j*ncells:(j+1)*ncells-1]  
    endfor
    ;resultant acceleration in each cell:
    agrav[*,xc,yc,zc]=total(acomponents,2)
    if (i mod 10000 eq 0) then print, i,'/',ncells
 endfor
;potentials and accelerations in cgs units:
cgsfactor=-G*mperH2*cellsize^2 ; for potential
Vgrav=cgsfactor*Vgrav
cgsfactor=G*mperH2*cellsize ; for acceleration
agrav=cgsfactor*agrav
; write binary cubes
OPENW,1,'hmm1_Vgrav.cube'
WRITEU,1,Vgrav
CLOSE,1
OPENW,2,'hmm1_agrav.cube'
WRITEU,2,agrav
CLOSE,2

end       
