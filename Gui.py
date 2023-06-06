#Gestione GUI
import tkinter as tk
from tkinter import filedialog,PhotoImage,Canvas
from PIL import Image, ImageTk
import os
import sys
import locale


import time
import threading
import logging
try:
    import tkinter as tk # Python 3.x
    import tkinter.scrolledtext as ScrolledText
except ImportError:
    import Tkinter as tk # Python 2.x
    import ScrolledText



def choose_directory():
    directory_path = filedialog.askdirectory()
    if directory_path:
        entry_var.set(directory_path)

def save_info():
    global datasetG2
    global SogliaCropG2
    global CHaddG2
    global FinalClosingG2
    global ProtrusG2
    global EdgesG2
    datasetG2 = dataset_entry.get()
    SogliaCropG2 = int(sogliaCrop_entry.get()) #Soglia per il crop del volume iniziale
    CHaddG2 = int(CHadd_entry.get())
    FinalClosingG2 = int(FinalClosing_entry.get())
    ProtrusG2=int(Protrus_entry.get())#Con ampiezza inferiore (3 o 4) peggiora
    EdgesG2=int(Edges_entry.get()) #(4 forse e' troppo)

    
    root.destroy()


def Gui(dataset, SogliaCrop, CHadd, FinalClosing, Protrus, Edges):
	
	global root
	root = tk.Tk()
	root.configure(bg="#33383B")
	root.title("KneeBones3Dify - DICOM SEGMENTATION")

	# Constructing the first frame, frame1
	frame1 = tk.LabelFrame(root, bg="#33383B",
	                    fg="#33383B", bd=0,pady=5)
	frame1.grid(row=0, column=0)
	frame2 = tk.LabelFrame(root, bg="#33383B",
	                    fg="#33383B",bd=0)
	frame2.grid(row=1, column=0)
	frame3 = tk.LabelFrame(root, bg="#33383B",
	                    fg="#33383B",bd=0)
	frame3.grid(row=2, column=0)

	media_img = Image.open("KneeBones3Dify_logo.png")
	media_img = media_img.resize((int(media_img.size[0]/4), int(media_img.size[1]/4)), Image.LANCZOS)
	media_imgI = media_img.resize((int(media_img.size[0]/3), int(media_img.size[1]/3)), Image.LANCZOS)
	media_img = ImageTk.PhotoImage(media_img)
	media_imgI = ImageTk.PhotoImage(media_imgI)
	root.wm_iconphoto(True, media_imgI)
	media_imgL = tk.Label(frame1, image=media_img, bg="#33383B",border=0,highlightthickness=0,borderwidth=0)
	media_imgL.grid(row=0,column=0)

	global entry_var
	entry_var = tk.StringVar(frame2,value=dataset)

	global dataset_entry
	dataset_entry = tk.Entry(frame2, textvariable=entry_var)
	dataset_entry.grid(row=1, column=2,sticky = "sw", padx=5)
	buttonData = tk.Button(frame2, text="DICOM Directory", command=choose_directory,fg = "#F8F9B5", bg="#33383B",anchor="e",justify="right")
	buttonData.grid(row=1, column=1,sticky = "sw", padx=5)

	global sogliaCrop_entry
	sogliaCrop_l = tk.Label(frame2, text="Intensity Threshold:",fg = "#F8F9B5", bg="#33383B",anchor="e",justify="right")
	sogliaCrop_l.grid(row=2, column=1,sticky = "sw", padx=5)
	sogliaCrop_entry = tk.Entry(frame2)
	sogliaCrop_entry.insert(600, str(SogliaCrop))
	sogliaCrop_entry.grid(row=2, column=2,sticky = "sw", padx=5)

	global CHadd_entry
	CHadd_l = tk.Label(frame2, text="Convex Hull Dilation:",fg = "#F8F9B5", bg="#33383B",anchor="e",justify="right")
	CHadd_l.grid(row=3, column=1, sticky = "sw", padx=5)
	CHadd_entry = tk.Entry(frame2)
	CHadd_entry.insert(6, str(CHadd))
	CHadd_entry.grid(row=3, column=2,sticky = "sw", padx=5)
	CHadd_l = tk.Label(frame2, text="",font=("Arial", 4),fg = "#F8F9B5", bg="#33383B")
	CHadd_l.grid(row=4, column=1)

	global FinalClosing_entry
	FinalClosing_l = tk.Label(frame2, text="Final Closing:", fg = "#F8F9B5", bg="#33383B",anchor="e",justify="right")
	FinalClosing_l.grid(row=1, column=3,sticky = "sw",padx=5)
	slider_var1 = tk.IntVar()
	slider_var1.set(FinalClosing)
	FinalClosing_entry = tk.Scale(frame2,fg = "#F8F9B5", bg="#33383B" ,from_=1, to=20, resolution=1,
	                    border=0,highlightthickness=0,borderwidth=0, 
	                    variable=slider_var1, orient=tk.HORIZONTAL)
	FinalClosing_entry.grid(row=1, column=4,sticky = "w",padx=5)

	global Protrus_entry
	Protrus_l = tk.Label(frame2, text="Protrusion Removal:",fg = "#F8F9B5", bg="#33383B",anchor="e",justify="right")
	Protrus_l.grid(row=2, column=3,sticky = "sw",padx=5)
	slider_var1 = tk.IntVar()
	slider_var1.set(Protrus)
	Protrus_entry =tk.Scale(frame2,fg = "#F8F9B5", bg="#33383B" ,from_=1, to=10, resolution=1,
	                    border=0,highlightthickness=0,borderwidth=0, 
	                    variable=slider_var1, orient=tk.HORIZONTAL)

	Protrus_entry.grid(row=2, column=4,sticky = "w",padx=5)

	global Edges_entry
	Edges_l = tk.Label(frame2, text="Final Dilation:",fg = "#F8F9B5", bg="#33383B",anchor="e",justify="right")
	Edges_l.grid(row=3, column=3,sticky = "sw",padx=5)
	slider_var = tk.IntVar()
	slider_var.set(Edges)
	Edges_entry =tk.Scale(frame2,fg = "#F8F9B5", bg="#33383B" ,from_=1, to=10, resolution=1,
	                    border=0,highlightthickness=0,borderwidth=0, 
	                    variable=slider_var, orient=tk.HORIZONTAL)
	Edges_entry.grid(row=3, column=4,sticky = "w",padx=5)


	submit_button = tk.Button(frame3, text="Ok", bg = "#1ED2F5", command=save_info)
	submit_button.grid(row=0, column=2, padx=15, pady=15)

	exit_button = tk.Button(frame3, text="Exit", bg = "#FB9764", command=exit)
	exit_button.grid(row=0, column=0, padx=15, pady=15)

	root.bind("<Return>", (lambda event: save_info()))
	root.resizable(False, False)
	root.mainloop()

	return datasetG2, SogliaCropG2, CHaddG2, FinalClosingG2, ProtrusG2, EdgesG2

def save_info2():
    global ris
    ris=True
    rootR.destroy()


def GuiFin():
	
	global rootR
	rootR = tk.Tk()
	rootR.configure(bg="#33383B")
	rootR.title("KneeBones3Dify - DICOM SEGMENTATION")

	# Constructing the first frame, frame1
	frame1 = tk.LabelFrame(rootR, bg="#33383B",
	                    fg="#33383B", bd=0,pady=5)
	frame1.grid(row=0, column=0)
	frame2 = tk.LabelFrame(rootR, bg="#33383B",
	                    fg="#33383B",bd=0)
	frame2.grid(row=1, column=0)
	frame3 = tk.LabelFrame(rootR, bg="#33383B",
	                    fg="#33383B",bd=0)
	frame3.grid(row=2, column=0)

	media_img = Image.open("KneeBones3Dify_logo.png")
	media_img = media_img.resize((int(media_img.size[0]/6), int(media_img.size[1]/6)), Image.LANCZOS)
	media_imgI = media_img.resize((int(media_img.size[0]/3), int(media_img.size[1]/3)), Image.LANCZOS)
	media_img = ImageTk.PhotoImage(media_img)
	media_imgI = ImageTk.PhotoImage(media_imgI)
	rootR.wm_iconphoto(True, media_imgI)
	media_imgL = tk.Label(frame1, image=media_img, bg="#33383B",border=0,highlightthickness=0,borderwidth=0)
	media_imgL.grid(row=0,column=0)


	global ris_entry
	ris_l = tk.Label(frame2, text="Would you like to continue with another segmentation?",fg = "#F8F9B5", bg="#33383B")
	ris_l.grid(row=3, column=1, padx=5, pady=5)

	submit_button = tk.Button(frame3, text="Yes", bg = "#1ED2F5", command=save_info2)
	submit_button.grid(row=0, column=3, padx=15, pady=15)
	exit_button = tk.Button(frame3, text="Exit", bg = "#FB9764", command=exit)
	exit_button.grid(row=0, column=2, padx=15, pady=15)

	rootR.bind("<Return>", (lambda event: save_info2()))
	rootR.resizable(False, False)
	rootR.mainloop()

	return ris

def error():
    rootE.destroy()

def GuiError():
	global rootE
	rootE = tk.Tk()
	rootE.configure(bg="#33383B")
	rootE.title("ERROR")

	# Constructing the first frame, frame1
	frame1 = tk.LabelFrame(rootE, bg="#33383B",
	                    fg="#33383B", bd=0,pady=5)
	frame1.grid(row=0, column=0)
	frame2 = tk.LabelFrame(rootE, bg="#33383B",
	                    fg="#33383B",bd=0)
	frame2.grid(row=1, column=0)
	frame3 = tk.LabelFrame(rootE, bg="#33383B",
	                    fg="#33383B",bd=0)
	frame3.grid(row=2, column=0)

	media_img = Image.open("error.png")
	media_img = media_img.resize((int(media_img.size[0]/6), int(media_img.size[1]/6)), Image.LANCZOS)
	media_imgI = media_img.resize((int(media_img.size[0]/3), int(media_img.size[1]/3)), Image.LANCZOS)
	media_img = ImageTk.PhotoImage(media_img)
	media_imgI = ImageTk.PhotoImage(media_imgI)
	rootE.wm_iconphoto(True, media_imgI)
	media_imgL = tk.Label(frame1, image=media_img, bg="#33383B",border=0,highlightthickness=0,borderwidth=0)
	media_imgL.grid(row=0,column=0)

	ris_l = tk.Label(frame2, text="File names information is empty. Cannot read series.\nPlease try again with a DICOM directory.",fg = "#F8F9B5", bg="#33383B")
	ris_l.grid(row=3, column=1, padx=5, pady=5)

	exit_button = tk.Button(frame3, text="Exit", bg = "#FB9764",command=error)
	exit_button.grid(row=0, column=2, padx=15, pady=15)

	rootE.bind("<Return>", (lambda event: error()))
	rootE.resizable(False, False)
	rootE.mainloop()
