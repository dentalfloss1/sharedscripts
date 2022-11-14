# I used this file to draw lines and shapes on top of an image 
# see astropyvisual.png for what this outputs


from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from scipy.special import binom
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True) # Enable LaTeX in matplotlib text
# Function below can be used either for the example data in the servers for astropy or for local files. 
filename = get_pkg_data_filename('/Users/schastain/Software/miscscripts/1582287955_sdp_l0.GRB200219A_im_3.fits') #'/media/sarah/Elements/GRB200219A/1582287955_sdp_l0.GRB200219A_im_3.fits') #'E:\\GRB200219A\\1582287955_sdp_l0.GRB200219A_im_3.fits')
hdu = fits.open(filename)[0] # Open the primary data, which ought to be an image
wcs = WCS(hdu.header, naxis=2) # parse header wcs







# RA,DEC, FOV radius. Specific values are from guessing and checking. 
uniquepoint = np.array([[343.0007591, -59.21352863, 0.25],[342.530408, -59.21438179, 0.25], [(1+0.0005)*342.7766486, -59.0499097, 0.25]],dtype=np.float64)
# uniquepoint = np.array([[1+90 , -1 , 1.5],[1-1.1+90, -1, 1.5], [(1-.55)+90, -(1.2)*1-.5, 1.5]],dtype=np.float64)

uniquesky = SkyCoord(ra=uniquepoint[:,0],dec=uniquepoint[:,1], unit='deg', frame='fk5') # Turn it into a SkyCoord obj
numrgns = len(uniquepoint)
# Line below creates an array of max possible number of regions. Don't really need to do this for just a figure. This is just because it is 
# taken directly from the simulations.
regions = np.zeros(np.uint32(numrgns + binom(numrgns,2) + binom(numrgns,3)), dtype={'names': ('ra', 'dec','identity', 'area', 'timespan', 'stop', 'start'), 'formats': ('f8','f8','U32','f8', 'f8', 'f8', 'f8')})
for i in range(numrgns): # Label the individual pointings. These regions are for example region 1 NOT 2 and 2 NOT 1 
    regions['identity'][i] = str(i)
    regions['ra'][i] = uniquepoint[i,0]
    regions['dec'][i] = uniquepoint[i,1]
    regions['area'][i] = (4*np.pi*np.sin(uniquepoint[i,2]*(np.pi/180/2))**2*(180/np.pi)**2) # Assumes single circular regions, for multiple pointings or other shapes this needs altering
    leftoff = i + 1
### SIMULATIONS HAVE COMMENTED OUT PART BECAUSE THERE WE ARE ABOUT START AND END TIMES OF OBSERVATIONS ###
# obssubsection = []
# for p in uniquepoint:
#     timeind = np.array([np.amax(np.argwhere((uniquepoint[:,0] == p[0]) & (uniquepoint[:,1] ==p[1]))), np.amin(np.argwhere((uniquepoint[:,0] == p[0]) & (uniquepoint[:,1] ==p[1])))])
#     matchregion = (regions['ra']==p[0]) & (regions['dec']==p[1]) & (regions['area']==(4*np.pi*np.sin(p[2]*(np.pi/180/2))**2*(180/np.pi)**2))
#     regions['timespan'][matchregion] = (observations['start'][timeind[0]] + observations['duration'][timeind[0]] - observations['start'][timeind[1]])
#     regions['stop'][matchregion] = observations['start'][timeind[0]] + observations['duration'][timeind[0]]
#     regions['start'][matchregion] = observations['start'][timeind[1]]
#     obssubsection.append([timeind[1],timeind[0],regions['identity'][matchregion][0]])
gamma = np.zeros((len(uniquepoint),len(uniquepoint)))
for i in range(len(uniquesky)-1): # Label intersections: For example: 1 AND 2
    for j in range(i+1,len(uniquesky)):
        if uniquesky[i].separation(uniquesky[j]).deg <= (uniquepoint[i,2] + uniquepoint[j,2]):
            d = uniquesky[i].separation(uniquesky[j]).rad
            r1 = uniquepoint[i,2]*np.pi/180.
            r2 = uniquepoint[j,2]*np.pi/180.
            gamma[i,j] = np.arctan((np.cos(r2)/np.cos(r1)/np.sin(d)) - (1/np.tan(d)))
            pa = uniquesky[i].position_angle(uniquesky[j])
            # https://arxiv.org/ftp/arxiv/papers/1205/1205.1396.pdf
            # and https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
            fullcone1 = 4*np.pi*np.sin(r1/2)**2 
            cutchord1 = 2*(np.arccos(np.sin(gamma[i,j])/np.sin(r1)) - np.cos(r1)*np.arccos(np.tan(gamma[i,j])/np.tan(r1))) 
            fullcone2 = 4*np.pi*np.sin(r2/2)**2
            cutchord2 = 2*(np.arccos(np.sin(gamma[i,j])/np.sin(r2)) - np.cos(r2)*np.arccos(np.tan(gamma[i,j])/np.tan(r2))) 
            centerreg = uniquesky[i].directional_offset_by(pa, gamma[i,j]*u.radian)
            regions['identity'][leftoff] = str(i)+'&'+str(j)
            regions['ra'][leftoff] = centerreg.ra.deg
            regions['dec'][leftoff] = centerreg.dec.deg
            regions['area'][leftoff] = (cutchord1 + cutchord2)*(180/np.pi)**2
            regions['start'][leftoff] = min(regions['start'][i],regions['start'][j])
            regions['stop'][leftoff] = max(regions['stop'][i],regions['stop'][j])
            regions['timespan'][leftoff] = regions['stop'][leftoff] - regions['start'][leftoff]
            leftoff+=1
scatterpointsra = []
scatterpointsdec = []
twonode = []
pointsorder = []
for i in range(len(uniquesky)-2): # repeat the above, but this time for triple overlapping regions
    for j in range(i+1,len(uniquesky)-1):
        for index3 in range(j+1,len(uniquesky)):
            if ((uniquesky[i].separation(uniquesky[j]).deg <= (uniquepoint[i,2] + uniquepoint[j,2])) and 
                (uniquesky[j].separation(uniquesky[index3]).deg <= (uniquepoint[j,2] + uniquepoint[index3,2])) and 
                (uniquesky[i].separation(uniquesky[index3]).deg <= (uniquepoint[i,2] + uniquepoint[index3,2]))):
                    r1 = uniquepoint[i,2]*np.pi/180
                    r2 = uniquepoint[j,2]*np.pi/180
                    r3 = uniquepoint[index3,2]*np.pi/180
                    # Get coordinates of the encircled(?) spherical triangle
                    # from the triangle formed between pointing center, overlap center, and overlap nodal point.

                    angle_offset = 90*u.deg
                    halfheightr4 = np.arccos(np.cos(r2)/np.cos(gamma[i][j])) 
                    point4key = np.where(regions['identity'] == str(i)+'&'+str(j))
                    point4ra = regions['ra'][point4key]
                    point4dec = regions['dec'][point4key]
                    point4sc = SkyCoord(ra=point4ra, dec=point4dec, unit='deg',frame='fk5')
                    point4pa = point4sc.position_angle(uniquesky[j])
                    point7sc = point4sc.directional_offset_by(point4pa + angle_offset, halfheightr4)
                    # print(point7sc.separation(uniquesky[index3]).deg)
                    if point7sc.separation(uniquesky[index3]).deg > uniquepoint[index3,2]:
                        point7sc = point4sc.directional_offset_by(point4pa - angle_offset, halfheightr4)

                    halfheightr5 = np.arccos(np.cos(r3)/np.cos(gamma[j][index3]))
                    point5key = np.where(regions['identity'] == str(j)+'&'+str(index3))
                    point5ra = regions['ra'][point5key]
                    point5dec = regions['dec'][point5key]
                    point5sc = SkyCoord(ra=point5ra, dec=point5dec, unit='deg',frame='fk5')
                    point5pa = point5sc.position_angle(uniquesky[index3])
                    point8sc = point5sc.directional_offset_by(point5pa + angle_offset, halfheightr5)
                    # print(point8sc.separation(uniquesky[i]).deg)
                    if point8sc.separation(uniquesky[i]).deg > uniquepoint[i,2]:
                        point8sc = point5sc.directional_offset_by(point5pa - angle_offset, halfheightr5)

                    halfheightr6 = np.arccos(np.cos(r1)/np.cos(gamma[i][index3]))
                    point6key = np.where(regions['identity'] == str(i)+'&'+str(index3))
                    point6ra = regions['ra'][point6key]
                    point6dec = regions['dec'][point6key]
                    point6sc = SkyCoord(ra=point6ra, dec=point6dec, unit='deg',frame='fk5')
                    point6pa = point6sc.position_angle(uniquesky[i])
                    point9sc = point6sc.directional_offset_by(point6pa + angle_offset, halfheightr6)
                    # print(point9sc.separation(uniquesky[j]).deg)
                    if point9sc.separation(uniquesky[j]).deg > uniquepoint[j,2]:
                        point9sc = point6sc.directional_offset_by(point6pa - angle_offset, halfheightr6)
                    scatterpointsra.extend([point7sc.ra,point8sc.ra,point9sc.ra])
                    scatterpointsdec.extend([point7sc.dec,point8sc.dec,point9sc.dec])
                    twonode.append([point6sc.directional_offset_by(point6pa + angle_offset, halfheightr6),point6sc.directional_offset_by(point6pa - angle_offset, halfheightr6)])
                    pointsorder.append([uniquesky[i],uniquesky[index3],uniquesky[j], point6sc[0], point9sc[0], point8sc[0], point7sc[0]]) # P1, P2, P3




# help(plt.subplot)
# wcs = WCS(None)
ax = plt.subplot(projection=wcs)
# remove the *0 to get the image in there. weird subscripting
# is necessary to make this work with some files
ax.imshow(hdu.data[0:][0:][0][0], vmin=-2.e-5, vmax=5.e-4, origin='lower',cmap='Greys')
ax.set_xlim(hdu.data.shape[3]/2 - 1200 , hdu.data.shape[3]/2 + 1000) # Set section of the image we want to see for x 
ax.set_ylim(hdu.data.shape[2]/2 - 1000, hdu.data.shape[2]/2 + 1000 ) # same but for y 
for p in uniquepoint: # plot our three spherical circles representing three pointings
    ax.add_patch(SphericalCircle((p[0] * u.deg, p[1] * u.deg), p[2] * u.degree,
                     edgecolor='red', facecolor='none',
                    transform=ax.get_transform('fk5'))) # ax.get_transform transforms to pixel coordinates
print("vertices: ", scatterpointsra, scatterpointsdec)
for t in twonode:
    ax.plot([t[0].ra, t[1].ra],[t[0].dec,t[1].dec], c='black',transform=ax.get_transform('fk5')) # Plot the line from P1 to P5
cnum = 1
# lines below make up the triangles. 
lines = np.array([[[pointsorder[0][0].ra.deg, pointsorder[0][1].ra.deg],[pointsorder[0][0].dec.deg,pointsorder[0][1].dec.deg]],
                  [[pointsorder[0][0].ra.deg, pointsorder[0][4].ra.deg], [pointsorder[0][0].dec.deg,pointsorder[0][4].dec.deg]],
                  [[pointsorder[0][1].ra.deg, pointsorder[0][4].ra.deg], [pointsorder[0][1].dec.deg,pointsorder[0][4].dec.deg]]])
# plot these lines 
for l in lines:
    ax.plot(l[0],l[1], c='black', transform=ax.get_transform('fk5'))

for p in pointsorder: # Need to label stuff
    for c in p:
        if cnum < 6:
            ax.scatter(c.ra,c.dec, s=10, marker='X', c='blue',transform=ax.get_transform('fk5'))
        annotationLocation = np.array(wcs.world_to_pixel(c)) # ax.annotate doesn't work with ax.get_transform, so we get pixel coords first 
        if cnum==1: # Do stuff at P1
            Aloc = np.copy(annotationLocation)
            Aloc += np.array([20,20]) # code like this just adjusts label positions
            annotationLocation += np.array([5,-70])
            ax.annotate("A",Aloc, c='black') # Annotate Angle A
        elif cnum==2:
            annotationLocation += np.array([-10,30])
            quadgoalpoint = np.copy(c) 
        elif cnum==3:
            annotationLocation += np.array([35,0])
        elif cnum==4:
            quadcorner = np.copy(c)
            Bloc = np.copy(annotationLocation)
            Bloc += np.array([10,-80])
            ax.annotate("B",Bloc, c='black') # annotate angle B
            annotationLocation += np.array([-100,40]) 
        elif cnum==5:
            Cloc = np.copy(annotationLocation)
            Cloc += np.array([-250,-30])
            ax.annotate("C",Cloc, c='black') # annotate angle C 
            annotationLocation += np.array([25,-5])
        if cnum < 6:
            ax.annotate("P"+str(cnum),annotationLocation, c='black') # annotate P1, P2, etc
        cnum = cnum + 1 

# We want to make a square to indicate a right angle
# One vertex is at P4, so we start there 
quadtilt = quadcorner.position_angle(quadgoalpoint) # and get the position angle to P2. 
quadp2 = quadcorner.directional_offset_by(quadtilt, 0.03*u.deg) # get second vertex by going up the line bc 0.03 deg
quadp3 = quadp2.directional_offset_by(quadtilt - 90*u.deg, 0.03*u.deg) # 90 degree turn and do it again 
quadp4 = quadp3.directional_offset_by(quadtilt + 180*u.deg, 0.03*u.deg) # 90 degree turn gives -180 now, do it again to close the box
ax.plot([quadp2.ra.deg, quadp3.ra.deg],[quadp2.dec.deg,quadp3.dec.deg], c='black',transform=ax.get_transform('fk5')) # plot one line of box
ax.plot([quadp3.ra.deg, quadp4.ra.deg],[quadp3.dec.deg,quadp4.dec.deg], c='black',transform=ax.get_transform('fk5')) # final line of box

# We now want to annotate the line labels. We do it by going up halfway up the line and making small adjustments to the label positions
for p in pointsorder:
    theta1sc = p[0].directional_offset_by(p[0].position_angle(p[4]),p[0].separation(p[4])/2.)
    theta1points = np.array(wcs.world_to_pixel(theta1sc))
    theta1points += np.array([-40,-80])
    ax.annotate("$\\theta_1$", theta1points, c='black')
    theta2sc = p[1].directional_offset_by(p[1].position_angle(p[4]),p[1].separation(p[4])/2.)
    theta2points = np.array(wcs.world_to_pixel(theta2sc))
    theta2points += np.array([0,5])
    ax.annotate("$\\theta_2$", theta2points, c='black')
    gamma1sc = p[0].directional_offset_by(p[0].position_angle(p[3]),p[0].separation(p[3])/2.)
    gamma1points = np.array(wcs.world_to_pixel(gamma1sc))
    gamma1points += np.array([-100,0])
    ax.annotate("$\\gamma_1$", gamma1points, c='black')
    gamma2sc = p[1].directional_offset_by(p[1].position_angle(p[3]),p[1].separation(p[3])/2.)
    gamma2points = np.array(wcs.world_to_pixel(gamma2sc))
    gamma2points += np.array([-90,20])
    ax.annotate("$\\gamma_2$", gamma2points, c='black')
    alphasc = p[3].directional_offset_by(p[3].position_angle(p[4]),p[3].separation(p[4])/2.)
    alphapoints = np.array(wcs.world_to_pixel(alphasc))
    alphapoints += np.array([-20,15])
    ax.annotate("$\\alpha$", alphapoints, c='black')

plt.xlabel("RA (J2000)")
plt.ylabel("Dec (J2000)")
plt.show()


plt.close()

ax2 = plt.subplot(projection=wcs)
ax2.imshow(hdu.data[0:][0:][0][0], vmin=-2.e-5, vmax=5.e-4, origin='lower',cmap='Greys')
ax2.set_xlim(hdu.data.shape[3]/2 - 1200 , hdu.data.shape[3]/2 + 1000) # Set section of the image we want to see for x 
ax2.set_ylim(hdu.data.shape[2]/2 - 1000, hdu.data.shape[2]/2 + 1000 ) # same but for y 
for p in uniquepoint: # plot our three spherical circles representing three pointings
    ax2.add_patch(SphericalCircle((p[0] * u.deg, p[1] * u.deg), p[2] * u.degree,
                     edgecolor='red', facecolor='none',
                    transform=ax.get_transform('fk5'))) # ax.get_transform transforms to pixel coordinates
cnum = 1

##### Change axes labels to RA/DEC  fk5 or whatever




for p in pointsorder: # Need to label stuff
    for c in p:
        if cnum >= 5:
            annotationLocation = np.array(wcs.world_to_pixel(c))
            if cnum==5:
                annotationLocation += np.array([20,-20])
            if cnum==6:
                annotationLocation += np.array([-70,-80])
            if cnum==7:
                annotationLocation += np.array([-10,20])
            ax2.annotate("P"+str(cnum), annotationLocation, c='black')
            ax2.scatter(c.ra,c.dec, s=10, marker='X', c='blue',transform=ax2.get_transform('fk5'))
            
        cnum += 1

lines2 = np.array([[[pointsorder[0][4].ra.deg, pointsorder[0][5].ra.deg],[pointsorder[0][4].dec.deg,pointsorder[0][5].dec.deg]],
                  [[pointsorder[0][5].ra.deg, pointsorder[0][6].ra.deg], [pointsorder[0][5].dec.deg,pointsorder[0][6].dec.deg]],
                  [[pointsorder[0][4].ra.deg, pointsorder[0][6].ra.deg], [pointsorder[0][4].dec.deg,pointsorder[0][6].dec.deg]]])
for l in lines2:
    ax2.plot(l[0],l[1], c='black', transform=ax.get_transform('fk5'))
plt.xlabel("RA (J2000)")
plt.ylabel("Dec (J2000)")
plt.show()
plt.close()
