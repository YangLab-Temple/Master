# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 12:40:19 2017

@author: Andrew Ruba


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import tkinter as tk
import numpy as np
import cv2
import os
import shutil
import tkinter.filedialog as tkfd

class splreader(tk.Frame):
    
    def __init__(self, master = None):
        tk.Frame.__init__(self, master)
        self.message = tk.StringVar()
        self.folder = None
        self.files = None
        self.mergedfiles = None
        self.pack()
        self.createwidgets()
        self.master = master
        
    def createwidgets(self):
        self.TOPCONTAINER = tk.Frame(self)
        self.TOPCONTAINER.pack({'side': 'top'})           
        
        self.LEFTTOPCONTAINER = tk.Frame(self.TOPCONTAINER)
        self.LEFTTOPCONTAINER.pack({'side': 'left'})   
        
        self.BOTTOMCONTAINER = tk.Frame(self)
        self.BOTTOMCONTAINER.pack({'side': 'bottom'})

        self.SELECT = tk.Button(self.LEFTTOPCONTAINER)
        self.SELECT['text'] = 'Select folder containing .spl files'
        self.SELECT['command'] = lambda: self.selectfolder()
        self.SELECT.pack({'side': 'top'})

        self.QUIT = tk.Button(self.LEFTTOPCONTAINER)
        self.QUIT['text'] = 'Quit'
        self.QUIT['fg'] = red
        self.QUIT['command'] = lambda: self.quit()
        self.QUIT.pack({'side': 'top'})

        self.MESSAGE = tk.Label(self.BOTTOMCONTAINER)
        self.MESSAGE['textvariable'] = self.message
        self.message.set('Select a folder.')
        self.MESSAGE.pack({'side': 'right'})

    def spl_to_tif(self, filename):
    
        ##converting .spl file to .tiff
        f = open(filename, 'rb')
    
        filetitle = filename.split('.')
        filetitle.pop()
        filetitle = '.'.join(filetitle)
        if os.path.exists(filetitle):
            shutil.rmtree(filetitle)
            os.mkdir(filetitle)
        else:
            os.mkdir(filetitle)
        filetitlename = filetitle.split('\\').pop()
    
        frame_header = ['f6','01','01','01','49','49','00','01']
        pic_width = 0
        pic_height = 0
        dimensions = []
        temp_header = []
        temp_frame = []
        frame_number = 1
        end_of_video = False
        looking_for_header = False
        found_header = False
        getting_pixel_data = False
        counter = 0
    
        while not end_of_video:
        
            try:
                val = str("{:02x}".format(ord(f.read(1))))
            except:
                end_of_video = True
                break
        
            if counter == 0 and val == frame_header[0] and not looking_for_header and not found_header and not getting_pixel_data:
                temp_header.append(val)
                getting_pixel_data = False
                looking_for_header = True
                continue

            if looking_for_header:
            
                could_be_header = True
                temp_header.append(val)
            
                for i in range(len(temp_header)):
                    if temp_header[i] != frame_header[i]:
                        could_be_header = False
                        looking_for_header = False
                        temp_header = []

                if could_be_header and len(temp_header) == len(frame_header):
                    found_header = True
                    looking_for_header = False
                    temp_header = []
                    counter = 248
                continue
                
            if found_header:
            
                if counter >= 235 and counter <= 238:
                    dimensions.append(val)
                    
                if counter == 1:
                    found_header = False
                    getting_pixel_data = True
                    pic_width = int(dimensions[1] + dimensions[0], 16)
                    pic_height = int(dimensions[3] + dimensions[2], 16)
                    dimensions = []
                    counter = pic_width * pic_height * 2
                
                counter -= 1
                continue

            if getting_pixel_data:

                temp_frame.append(val)
            
                if counter == 0:
                    self.message.set(str(frame_number) + ' frames converted from ' + filetitlename + '.spl')
                    self.update()
                    getting_pixel_data = False
                    temp_frame = [ x+y for x,y in zip(temp_frame[1::2], temp_frame[0::2]) ]
                    temp_frame = [ int(x, 16) for x in temp_frame ]
                    temp_frame = np.reshape(temp_frame, (pic_height, pic_width))
                    temp_frame = np.array(temp_frame).astype(np.uint16)
                    cv2.imwrite(filetitle + '/' + str(frame_number) + '.tif', temp_frame)
                    frame_number += 1
                    temp_frame = []
                    continue
        
                counter -= 1
                continue

    def mergefiles(self):

        if self.files != None:

            self.SELECT['text'] = 'Select folder containing .spl files'
            self.SELECT['fg'] = black
            self.SELECT['command'] = lambda: self.selectfolder()   
                
            for filename in self.files:
                self.spl_to_tif(filename)
            self.message.set('Done. Select a folder.')

    def listfiles(self, folder, file_ext):

        folder_list = os.listdir(folder)

        file_list = []
    
        for file_path in folder_list:
            split_ext = file_path.split(os.extsep).pop()
            if file_ext == split_ext:
                file_list.append(os.path.normpath(os.path.join(folder, file_path)))
    
        return file_list

    def selectfolder(self):
        self.folder = tkfd.askdirectory(initialdir = 'C:\\', parent = self, title = 'Select a folder', mustexist = True)
        self.files = self.listfiles(self.folder, 'spl')
        
        self.message.set(str(len(self.files)) + ' .spl file(s) found.')
        
        if len(self.files) > 0:
            self.SELECT['text'] = 'convert .spl files'
            self.SELECT['fg'] = green
            self.SELECT['command'] = lambda: self.mergefiles()
        
        print (self.files)

root = tk.Tk()
root.geometry('250x75')
root.title(string = '.spl reader')

orange = '#FF7F22'
black = '#000000'
gray = '#D3D3D3'
white = '#FFFFFF'
red = '#FF0000'
green = '#008000'

app = splreader(master = root)
app.mainloop()
root.destroy()