# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 22:58:09 2018

@author: Roman
"""

import numpy as np
import cv2
import math
import matplotlib.pyplot as plt
import tkinter as tk

#initialization of variables used for calculation of velocity and acceleration
count = 0
acelList = []
poldest=0
velList = []
oldDisp = []

# input video
cap = cv2.VideoCapture('cellmoves.mp4')

#cap = cv2.VideoCapture('cellmoves.mp4')
fps = cap.get(cv2.CAP_PROP_FPS)


# params for ShiTomasi corner detection
# maxCorners = maximum # of corners to return. Limits by strongest
#quality level= Minimum accepted quality level. Minimum eigenvalue
#minDistance = Minimum Euclidean distance between corners
#Blocksize = Avg. block for computing derivative covariation matrix over each pixel neighborhood
feature_params = dict( maxCorners = 100,
                       qualityLevel = 0.3,
                       minDistance = 7,
                       blockSize = 7 )

# Parameters for lucas kanade optical flow
# winSize --> size of the search window at each pyramid level
# maxLevel --> for pyramidal implementation. If 0 , pyramids are not (one level). If 1, 2 levels to pyramid, etc.
#criteria--> termination criteria of iterative search algorithm
#Either when window moves by less than epsilon (10) or max number of iterations (0.03)
lk_params = dict( winSize  = (15,15),
                  maxLevel = 2,
                  criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 10, 0.03))

# Create some random colors
color = np.random.randint(0,255,(100,3))

# Get return value (0 or 1) and frame image
ret, old_frame = cap.read()

#only works if return value is true
if ret:
    #cvtColor converts from one color space to another
    #cv libraries use the BGR colour space, so need to convert RGB to BGR
    old_gray = cv2.cvtColor(old_frame, cv2.COLOR_BGR2GRAY)
    
    #find initial corner locations
    #mask is optional region of interest. May specify region where corners are detected
    #needs to be in CV_8UC1 type
    p0 = cv2.goodFeaturesToTrack(old_gray, mask = None, **feature_params)

    # Create a mask image for drawing purposes
    #mask is declared as array of zeros with same size as image array
    #this array will later be used to draw the paths behind the points
    mask = np.zeros_like(old_frame)

#arrays for max values (vel and/or acel) and time to be graphed
xaxis=[]
yaxis=[]
    
while(1):
    #arrays for calculating current change in displacement, velocity and accel. of all points
    curDisp = []
    curVels=[]
    curAcel =[]
    
    ret,frame = cap.read()
    
    if ret:
        frame_gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

    # calculate optical flow
    #p1 is output vector of points; 
    #st is status-->1 if flow is found, 0 if not. Read as output vector
    #err-> output vector of errors
        p1, st, err = cv2.calcOpticalFlowPyrLK(old_gray, frame_gray, p0, None, **lk_params)


    # Creates a cleaner array which only returns parts of p1 and p2 where st = 1
        good_new = p1[st==1]
        good_old = p0[st==1]
                
        #calculate euclidean distance for all good points, place into array of displacements
        for i in range(0, len(good_old)):
            xdif = abs(good_old[i][0]-good_new[i][0])
            ydif = abs(good_old[i][1]-good_new[i][1])
            disp = math.sqrt(abs(math.pow(xdif, 2) + math.pow(ydif, 2)))
            curDisp.append(disp)
        
        
    # draw the tracks for every point's movement
        for i,(new,old) in enumerate(zip(good_new,good_old)):
            a,b = new.ravel()
            c,d = old.ravel()
            mask = cv2.line(mask, (a,b),(c,d), color[i].tolist(), 2)
            frame = cv2.circle(frame,(a,b),5,color[i].tolist(),-1)
            
        #
        img = cv2.add(frame,mask)        
        
    #show current frame (loop works fast, making video )
        cv2.imshow('frame',img)
        #k = cv2.waitKey(30) & 0xff
        if cv2.waitKey(1)==13:
            break
        
        #to create a list for change in displacement. e.g. velocity
        if (count ==0):
            velList = curDisp
            
        else:
            for i in range(0, len(good_old)):
                if (curDisp[i]):
                    velList[i] += curDisp[i]
                else:
                    velList[i] = velList[i]        
        
    # Now update the previous frame and previous points
        old_gray = frame_gray.copy()
        p0 = good_new.reshape(-1,1,2)       
        
        #to calculate change in velocity(change in change in displacement, d^2x). e.g. acceleration
        if (count !=0):
            for i in range(0, len(good_old)):
                
                doubleDispDiff = abs(oldDisp[i]-curDisp[i])
                
                curAcel.append(doubleDispDiff)
                
                if (count==1):
                    acelList.append (doubleDispDiff)
                    
                else:
                    acelList[i] += doubleDispDiff
        
        poldest = good_old      
        oldDisp = curDisp
        
      #trying to calculate current maximum velocity and acceleration based on frame
        if (count>1):
            curVels = curDisp
            curAcel = acelList
            sortedCurV=sorted(curVels, reverse=True)
            sortedCurA=sorted(curAcel, reverse=True)
            
            fastCurV = sortedCurV[0]
            fastCurV *= float(fps)
            fastCurV= ("%.4f" % fastCurV)
            
            fastCurA = sortedCurA[0]
            curSec = count/fps
            fastCurA *= float(fps*fps)
            fastCurA= ("%.4f" % fastCurA)
           # print(fastCurA)
            
            curSec = ("%.2f" % curSec)
            
            #print ("at "+ str(curSec) + " seconds: ")
            #print("velocity: "+ str(fastCurV)+ "pix/sec.")
            #print("acceleration: " + str(fastCurA)+ " pix/sec^2")
            xaxis.append(curSec)
            yaxis.append(fastCurA)            
            
        count+=1       
                
    #if no more frames are detected, exit the loop  
    else:
        break

#destroy window and release video play
cv2.destroyAllWindows()
cap.release()

#find average velocity and acceleration through all frames 
for i in range (0, len(velList)):
    velList[i] = velList[i]/(count-1)
    
for i in range (0, len(acelList)):  
    acelList[i] = acelList[i]/(count-2)

#sort velocity and acceleration from biggest to smallest values
velList = sorted(velList, reverse=True)
acelList = sorted(acelList, reverse=True)

#we only care about the biggest value, which is the first one
fastestVel = ("%.2f" % (velList[0]))
fastestAcel = ("%.2f" % (acelList[0]))

#frame * frame/sec = sec
seconds = count/fps

print ("TOTAL VIDEO TIME: " + str(seconds))

print ("The average displacement (pix/frame) for the fastest object is: " + str(fastestVel) + "pixels/frame")
print ("The average difference in displacement (pix/frame) for the fastest object is: " + str(fastestAcel))

#velocity is currently as pixels/frame. This is multiplied by frames/sec to give us pixels/sec
realVel = ("%.2f" % (float(fastestVel) * fps))
print ("the real average velocity of the fastest object was: " + str(realVel) + "pixels/sec")

realAcel = ("%.2f" % (float(fastestAcel) * fps * fps))
print ("the real average acceleration of the fastest object was: " + str(realAcel) + "pixels/sec^2")

plt.title("Acceleration vs. Time")
plt.xlabel("Time")
plt.ylabel("Acceleration")
plt.plot(xaxis, yaxis, linewidth=2.0)


root = tk.Tk()


textLabel = tk.Label(root, text="The Average Velocity of the fastest object was: "+ str(realVel) + "pixels/sec")
textLabel.pack()

#root.mainloop()

#print (fps)
#print (seconds)

    