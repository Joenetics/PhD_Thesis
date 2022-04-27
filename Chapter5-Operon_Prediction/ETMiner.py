from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from guizero import App, PushButton, Text, TextBox, warn, ListBox, MenuBar, Window
import os
import json
from datetime import datetime
from random import randrange
import numpy as np
import math
import matplotlib.pyplot as plt
#import pkg_resources.py2_warn
import time
import pandas as pd
import glob
import openpyxl
from openpyxl.utils import get_column_letter
from openpyxl.styles import Color, PatternFill, Font, Border
from distutils.core import setup
from PIL import ImageTk  # package 'Pillow'
from PIL import Image
import tkinter
from tkinter import filedialog

global Sequences, Haem_Weight_Dict, operon_type_dict, print_dict, count_dict, top_hits, TM_thresolh, CxxH_threshold, \
    output_file_name, input_csv, TM_File, cxxh_file, plt, output_directory_images
global input_csv
global cxxh_file
global TM_File
global output_dir
global currdir
datenow = str(datetime.now()).split(" ")[0]

python_testing = False  # True == when using code with .py, False == using code as .exe
if python_testing:
    currdir = os.getcwd() + "\\"
    try:
        os.mkdir(os.getcwd() + "\\Output_Files\\")
    except:
        pass
    output_dir = currdir + "\\Output_Files\\"
    output_directory_images = output_dir + "Custom_Operon_images\\"
    sorted_variables_directory = currdir + "Stored_Variables\\"
else:
    currdir = '\\'.join(os.getcwd().split("\\")[:-2]) + "\\"
    try:
        os.mkdir('\\'.join(os.getcwd().split("\\")[:-2]) + "\\Output_Files\\")
    except:
        pass
    output_dir = currdir + "\\Output_Files\\"
    output_directory_images = output_dir + "Custom_Operon_images\\"
    sorted_variables_directory = currdir + "Stored_Variables\\"


#Electron Transfer (data)Miner
#ETMiner !
app = App(title="ETMiner: Electron Transfer Miner", layout="grid", height=320, width=700)  ## Title of GUI. Made it be the name of program
app.tk.iconbitmap("ETMiner.ico")
window = Window(app, title="Add custom operon searches (from CSV output)", layout="grid", height=150, width=350)
window2 = Window(app, title="Create operon sequence images", layout="grid", height=250, width=400)
window.hide()
window2.hide()
global customs
customs = []


def update_status(new_status, colour):
    newstring = "Status: %s" % str(new_status)  # put newline every 24 chars
    status.value = '-\n'.join(newstring[i:i + 20] for i in range(0, len(newstring), 20))
    if colour:
        status.bg = colour
    else:
        status.bg = "White"
    # print("Status: %s" % str(new_status), file = open("out.txt", "a+"))


def scatter_maker(output_location, tm_thresh, cxxh_thresh, scatter_data,  tm_thresh_upper, cxxh_thresh_upper):
    print(f"{datetime.now()}: Started making scatterplots.")
    #time.sleep(30)
    global plt
    output_location_scatter = "%sTM%s_%s_Haem%s_%s\\Scatterplots\\" % (
        output_location, str(tm_thresh), str(tm_thresh_upper), str(cxxh_thresh), str(cxxh_thresh_upper))
    try:
        os.mkdir(output_location_scatter)
    except:
        pass
    if scatter_data:
        major_csv = {}

        for operon_type in scatter_data:
            if "raws" in operon_type:
                all_x = []
                all_y = []
                all_names = []
                data_in_type = scatter_data[operon_type]
                data_in_names = scatter_data[operon_type.replace("raws", "operonname")]
                for datapoint in range(0, len(data_in_type)):
                    operon_results = data_in_type[datapoint]
                    for each_protein in operon_results:
                        data = each_protein.split(":")
                        if float(upper_limit.value) >= float(data[1]) >= float(lower_limit.value):
                            genename = ",".join(data_in_names[datapoint].split(",")[0:3])
                            if "" in genename:
                                all_names.append(genename.replace("", ""))
                            else:
                                all_names.append(genename)
                            all_x.append(int(data[0]))
                            all_y.append(float(data[1]))
                #print(f"*****\n{all_x}\n{all_y}")
                major_csv[operon_type.strip("_raws") + "_OperonNames"] = all_names
                del all_names
                major_csv[operon_type.strip("_raws") + "_X_Axis"] = all_x
                major_csv[operon_type.strip("_raws") + "_Y_Axis"] = all_y
                colors = (0, 0, 0)
                area = np.pi * 3
                # Plot
                all_y = [math.log10(i) for i in all_y]
                x = np.array(all_x)
                y = np.array(all_y)
                if not all_x:
                    all_x=[0]
                    all_y=[0]
                x_max = int(max(all_x))
                step = int(float(x_max)/10)
                y_max = int(max(all_y))
                step_y = int(float(y_max)/10)
                if step_y == 0:
                    step_y = 1
                if step == 0:
                    step = 1
                #plt.xticks(range(len(tickers)), tickers)
                plt.scatter(x, y, s=4) # (s = int) # for size of dots
                #plt.xticks(range(min(new_x), max(all_x)+1, step))
                plt.yticks(range(0, int(max(all_y)) + step_y, step_y))
                plt.xticks(range(int(min(all_x)), int(max(all_x)) + step, step))
                #print(f"For operon type {operon_type}, we get {list(range(min(new_x), max(all_x), step))}")
                #plt.scatter(x, y, s=area, c=colors, alpha=0.5)
                plt.title(f'Scatter plot for {operon_type}')
                plt.xlabel('Haem Number')
                plt.ylabel('Log10 Protein kDa/Haem')
                #plt.show()
                #print('%s\\%s.png' % (output_location_scatter, operon_type))
                plt.savefig('%s\\%s.png' % (output_location_scatter, operon_type), dpi=300, bbox_inches="tight")
                # plt.close()
                plt.clf()
                time.sleep(0.2)
                #print('%s\\%s.png' % (output_location_scatter, operon_type))
        longest_col = 0
        for col in major_csv:
            if len(major_csv[col]) > longest_col:
                longest_col = len(major_csv[col])
        for col in major_csv:
            if len(major_csv[col]) < longest_col:
                major_csv[col] = major_csv[col] + ([""]*(longest_col-len(major_csv[col])))

        output_df = pd.DataFrame(major_csv)

        output_df.to_excel(f"{output_location_scatter}\\All_Scatterplot_Data.xlsx", index=False, sheet_name="Data_For_Scatterplots")
        del output_df
        del major_csv
    print(f"{datetime.now()}: Finished making scatterplots.")
    #time.sleep(30)
    print(f"{datetime.now()}: All done.")
    app.after(1000, update_status, args=["Finished", "Green"])


def histo_maker(output_location, tm_thresh, cxxh_thresh, all_data,  tm_thresh_upper, cxxh_thresh_upper, scatterplot_data):
    print(f"{datetime.now()}: Started making histogram data.")
    #time.sleep(30)
    try:
        global plt
        output_location_hist = "%sTM%s_%s_Haem%s_%s\\Histograms\\" % (output_location, str(tm_thresh),
                                                                      str(tm_thresh_upper),
                                                                      str(cxxh_thresh), str(cxxh_thresh_upper))
        try:
            os.mkdir(output_location_hist)
        except:
            pass
        histogram_data_dict = {"Cyc2": [], "MtrAB": [], "MtrCAB": [], "custom": [], "Other": []}
        # ProteinID-Cxxh_motifs-TM_regions-MolWeight
        for line in all_data:
            data = line.split(",")
            present_types = data[0].split(";")
            #print(data)
            #sys.exit()
            haems_present = list()
            tms_present = list()
            mw_present = list()
            for prot in range(4, len(data)):
                haems_present.append(int(data[prot].split("-")[1]))
                tms_present.append(int(data[prot].split("-")[2]))
                mw_present.append(float(data[prot].split("-")[3]))
            for prot_type in present_types:
                if prot_type == "Cyc2":
                    for i in range(0, len(tms_present)):
                        if tms_present[i] >= tm_thresh and haems_present[i] >= cxxh_thresh and tms_present[i] <= tm_thresh_upper and haems_present[i] <= cxxh_thresh_upper:
                            histogram_data_dict["Cyc2"].append("%s,%s,%s" % (data[1], str(tms_present[i]), str(haems_present[i])))
                            #if "Cyc2" in histogram_data_dict:
                            #    histogram_data_dict["Cyc2"].append("%s,%s,%s" % (data[1], str(tms_present[i]), str(haems_present[i])))
                            #    break
                            #else:
                            #    histogram_data_dict["Cyc2"] = ["%s,%s,%s" % (data[1], str(tms_present[i]), str(haems_present[i]))]
                            #    break
                elif prot_type == "MtrAB":
                    for i in range(0, len(tms_present)):
                        if i < (len(tms_present)-1) and tms_present[i] >= tm_thresh and haems_present[i+1] >= cxxh_thresh and tms_present[i] <= tm_thresh_upper and haems_present[i+1] <= cxxh_thresh_upper:
                            histogram_data_dict["MtrAB"].append("%s,%s,%s" % (data[1], str(tms_present[i]), str(haems_present[i + 1])))
                            #if "MtrAB" in histogram_data_dict:
                            #    histogram_data_dict["MtrAB"].append("%s,%s,%s" % (data[1], str(tms_present[i]), str(haems_present[i+1])))
                            #else:
                            #    histogram_data_dict["MtrAB"] = ["%s,%s,%s" % (data[1], str(tms_present[i]), str(haems_present[i+1]))]
                            #break
                        elif i > 0 and tms_present[i] >= tm_thresh and haems_present[i-1] >= cxxh_thresh and tms_present[i] <= tm_thresh_upper and haems_present[i-1] <= cxxh_thresh_upper:
                            histogram_data_dict["MtrAB"].append("%s,%s,%s" % (data[1], str(tms_present[i]), str(haems_present[i - 1])))
                            #if "MtrAB" in histogram_data_dict:
                            #    histogram_data_dict["MtrAB"].append("%s,%s,%s" % (data[1], str(tms_present[i]), str(haems_present[i - 1])))
                            #else:
                            #    histogram_data_dict["MtrAB"] = ["%s,%s,%s" % (data[1], str(tms_present[i]), str(haems_present[i - 1]))]
                            #break
                elif prot_type == "MtrCAB":
                    for i in range(0, len(tms_present)):
                        if i < (len(tms_present)-1)  and i > 0 and tms_present[i] >= tm_thresh and haems_present[i+1] >= cxxh_thresh\
                                and haems_present[i-1] >= cxxh_thresh and tms_present[i] <= tm_thresh_upper and haems_present[i+1] <= cxxh_thresh_upper and haems_present[i-1] <= cxxh_thresh_upper:
                            histogram_data_dict["MtrCAB"].append("%s,%s,%s %s" % (data[1], str(tms_present[i]), str(haems_present[i + 1]), str(haems_present[i - 1])))
                            #if "MtrCAB" in histogram_data_dict:
                            #    histogram_data_dict["MtrCAB"].append("%s,%s,%s %s" % (data[1], str(tms_present[i]), str(haems_present[i+1]), str(haems_present[i-1])))
                            #else:
                            #    histogram_data_dict["MtrCAB"] = ["%s,%s,%s %s" % (data[1], str(tms_present[i]), str(haems_present[i+1]), str(haems_present[i-1]))]
                            #break
                        elif i > 1 and tms_present[i] >= tm_thresh and haems_present[i-1] >= cxxh_thresh and haems_present[i-2] >= cxxh_thresh and \
                                tms_present[i] <= tm_thresh_upper and haems_present[i-1] <= cxxh_thresh_upper and haems_present[i-2] <= cxxh_thresh_upper:
                            histogram_data_dict["MtrCAB"].append("%s,%s,%s %s" % (data[1], str(tms_present[i]), str(haems_present[i - 1]), str(haems_present[i - 2])))

                            #if "MtrCAB" in histogram_data_dict:
                            #    histogram_data_dict["MtrCAB"].append("%s,%s,%s %s" % (data[1], str(tms_present[i]), str(haems_present[i - 1]), str(haems_present[i - 2])))
                            #else:
                            #    histogram_data_dict["MtrCAB"] = [
                            #        "%s,%s,%s %s" % (data[1], str(tms_present[i]), str(haems_present[i - 1]), str(haems_present[i - 2]))]
                            #break
                        elif i < (len(tms_present)-2) and tms_present[i] >= tm_thresh and haems_present[i+1] >= cxxh_thresh and haems_present[i+2] >= cxxh_thresh and \
                                tms_present[i] <= tm_thresh_upper and haems_present[i+1] <= cxxh_thresh_upper and haems_present[i+2] <= cxxh_thresh_upper:
                            histogram_data_dict["MtrCAB"].append("%s,%s,%s %s" % (data[1], str(tms_present[i]), str(haems_present[i + 1]), str(haems_present[i + 2])))
                            #if "MtrCAB" in histogram_data_dict:
                            #    histogram_data_dict["MtrCAB"].append("%s,%s,%s %s" % (data[1], str(tms_present[i]), str(haems_present[i + 1]), str(haems_present[i + 2])))
                            #else:
                            #    histogram_data_dict["MtrCAB"] = [
                            #        "%s,%s,%s %s" % (data[1], str(tms_present[i]), str(haems_present[i + 1]), str(haems_present[i + 2]))]
                            #break
                else:
                    tm_bit = list()
                    hm_bit = list()
                    for i in range(0, len(tms_present)):
                        if tms_present[i] >= tm_thresh and tms_present[i] <= tm_thresh_upper:
                            tm_bit.append(str(tms_present[i]))
                        if haems_present[i] >= cxxh_thresh and haems_present[i] <= cxxh_thresh_upper:
                            hm_bit.append(str(haems_present[i]))
                    histogram_data_dict[prot_type].append("%s,%s,%s" % (data[1], " ".join(tm_bit), " ".join(hm_bit)))
                    #if prot_type in histogram_data_dict:
                    #    histogram_data_dict[prot_type].append("%s,%s,%s" % (data[1], " ".join(tm_bit), " ".join(hm_bit)))
                    #else:
                    #    histogram_data_dict[prot_type] = ["%s,%s,%s" % (data[1], " ".join(tm_bit), " ".join(hm_bit))]
        print(f"{datetime.now()}: Histogram data ready for heatmaps and histograms.")
        #time.sleep(30)
        to_del = []
        for type in histogram_data_dict:
            if not histogram_data_dict[type]:
                to_del.append(type)
        for ea in to_del:  # remove empty plots, to not mess up future data printing.
            del histogram_data_dict[ea]
        output_location_heat = "%sTM%s_%s_Haem%s_%s\\Heatmaps\\" % (
        output_location, str(tm_thresh), str(tm_thresh_upper), str(cxxh_thresh), str(cxxh_thresh_upper))
        try:
            os.mkdir(output_location_heat)
        except:
            pass
        for type in histogram_data_dict:
            temp_csv_name = "%sJoinedData%s_TMs_%s_%s_Haems_%s_%s.csv" % (output_location_hist, type, str(tm_thresh), str(tm_thresh_upper), str(cxxh_thresh), str(cxxh_thresh_upper))
            print("Taxon ID,TMs,Haems\n", "\n".join(histogram_data_dict[type]), sep="", file=open(temp_csv_name, "w"))
            #temp_csv = open("%s\\JoinedData%s_TMs_%s_%s_Haems_%s_%s.csv" % (output_location_hist, type, str(tm_thresh), str(tm_thresh_upper), str(cxxh_thresh), str(cxxh_thresh_upper)), "w")
            #print("Taxon ID,TMs,Haems", sep="", file=temp_csv)
            #print("\n".join(histogram_data_dict[type]), sep = "", file=temp_csv)
            #temp_csv.close()  # send data for manual print somewhere.
            heat_xlsx = {}
            histo_dict = {"All": []}
            heat_dict = dict()
            haems_in_histo = []
            lowest_value_haem = cxxh_thresh
            lowest_value_tm = tm_thresh
            highest_value_haem = cxxh_thresh
            highest_value_tm = tm_thresh
            for entry in histogram_data_dict[type]:
                tms = entry.split(",")[1].split(" ")
                haems = entry.split(",")[2].split(" ")
                for each_tm in tms:
                    if int(each_tm) > highest_value_tm:
                        highest_value_tm = int(each_tm)
                    for each_haem in haems:
                        haems_in_histo.append(int(each_haem))
                        if int(each_haem) > highest_value_haem:
                            highest_value_haem = int(each_haem)
                        if int(each_tm) in histo_dict:
                            histo_dict[int(each_tm)].append(int(each_haem))
                        elif haems != [""]:
                            histo_dict[int(each_tm)] = [int(each_haem)]
                        histo_dict["All"].append(int(each_haem))
                        #if "All" in histo_dict:
                        #    histo_dict["All"].append(int(each_haem))
                        #else:
                        #    histo_dict["All"] = [int(each_haem)]
            #full_range_haem = list(range(lowest_value_haem, highest_value_haem + 1))

            full_range_haem = list(sorted(set(haems_in_histo)))

            full_range_tm = list(range(lowest_value_tm, highest_value_tm + 1))

            for each_tm in full_range_tm:
                if each_tm in histo_dict:
                    #print("Type is %s and TM is %s" % (type, each_tm))
                    heat_dict[each_tm] = [0] * len(full_range_haem)
                    for each_haem in histo_dict[each_tm]:
                        # basically, find where haem fits into heat map, then add one to said point in heat map.
                        #heat_dict[each_tm][full_range_haem.index(each_haem)] = heat_dict[each_tm][full_range_haem.index(each_haem)] + 1
                        heat_dict[each_tm][full_range_haem.index(each_haem)] += 1
                #else:  # if TM doesnt exist in the data, just make it blank on heatmap. delete this 'else' to remove blanks
                #    heat_dict[each_tm] = [0] * len(full_range_haem)
            for list_key in heat_dict:
                #print("******")
                #print(heat_dict[list_key])
                for value in range(0, len(heat_dict[list_key])):
                    #print("Current index: %s" % value)
                    if heat_dict[list_key][value] == 0:
                        heat_dict[list_key][value] = 0
                    else:
                        heat_dict[list_key][value] = math.log10(heat_dict[list_key][value])
                #print(heat_dict[list_key])
            marks = np.array(list(reversed(list(heat_dict.values()))))
            heat_x = full_range_haem
            # heat_y = full_range_tm
            heat_y = list(reversed(list(heat_dict.keys())))
            #plt.xticks(heat_x)
            #plt.xticks(ticks=np.arange(len(heat_x)), labels=heat_x)
            plt.xticks(ticks=np.arange(len(heat_x)), labels=heat_x, rotation=90)
            plt.yticks(ticks=np.arange(len(heat_y)), labels=heat_y)
            hm = plt.imshow(marks, cmap=operon_choice.value, interpolation="nearest")  #cmap = blues,hot, cool,.. BLUES IS BEST
            plt.colorbar(hm)
            plt.title("Heatmap of %s Proteins" % type)
            plt.xlabel("Haem number in operon")
            plt.ylabel("TM number in operon")
            plt.savefig('%s\\%s.png' % (output_location_heat, type), dpi=300, bbox_inches="tight")
            #heat_xlsx[f"{type}_X_axis"] = heat_x
            #output_df = pd.DataFrame(major_csv)
            #output_df.to_excel(f"{output_location_scatter}\\All_Scatterplot_Data.xlsx", index=False, sheet_name="Data_For_Scatterplots")
            #make xlsx here4
            # plt.close()
            plt.clf()
            time.sleep(0.2)

            for each_plot in histo_dict:
                x_values = histo_dict[each_plot]
                x_axis = list(range(0, max(x_values)))
                step = int(float(max(x_axis))/10)
                if step == 0:
                    step = 1
                #plt.xticks(range(len(tickers)), tickers)
                """
                if each_plot == 48:
                    print("***")
                    print(x_values)
                    print(x_axis)
                    print(new_x)
                    print(range(min(new_x), max(x_values), step))"""
                try:
                    plt.hist(x_values, color=operon_choice.value.replace("s", ""), rwidth=0.95, bins=np.arange(max(x_axis))-0.5)

                except:  # if colour doesn't exist, use default (blue)
                    plt.hist(x_values)
                if step <= 1:
                    plt.xticks(np.arange(0, max(x_axis), 1))
                else:
                    plt.xticks(np.arange(0, max(x_axis), step))
                #plt.xticks(range(min(x_values), max(x_values) + step, step))
                #plt.xlim([x_axis[0]-1, max(x_values) + 1])
                plt.title("%s Proteins with %s TM regions" % (type, str(each_plot)))
                plt.xlabel("Haems in operon protein")
                plt.ylabel("Number of species")
                plt.savefig('%s%s_%sTMs.png' % (output_location_hist, type, str(each_plot)), dpi=300, bbox_inches="tight")
                #plt.close()
                plt.clf()
                time.sleep(0.2)

        plt.close('all')
        #app.after(500, update_status, args=["Finished", "Red"])
        #app.after(10000, redo_plt)
    except Exception as e:
        print(e)

    app.after(500, update_status, args=["Working on scatter plots...", "Yellow"])
    try:
        app.after(1300, scatter_maker, args=[output_dir, TM_thresolh, CxxH_threshold, scatterplot_data, TM_thresolh_upper, CxxH_threshold_upper])
    except:
        app.after(1000, update_status, args=["Failed to generate scatter plots", "red"])

def redo_plt():
    global plt
    import matplotlib.pyplot as plt

def fig_maker(file, filetype, image_dir):
    try:
        record = SeqIO.read("%s%s.gb" % (image_dir, file), "genbank")
        #record = SeqIO.read("NC_005816.gb", "genbank")
        #gd_diagram = GenomeDiagram.Diagram("Yersinia pestis biovar Microtus plasmid pPCP1")
        gd_diagram = GenomeDiagram.Diagram("Bacterial Operon on Shotgun Genome Sequence")
        gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
        gd_feature_set = gd_track_for_features.new_set()

        Highest_Number = 1
        Lowest_Number = 9999999999
        for feature in record.features:
            if feature.type != "gene":
                #Exclude this feature
                continue

            protein = (feature.qualifiers["locus_tag"][0]).split("-")
            if int(str(feature.location).split("]")[0].split(":")[1]) > Highest_Number:
                Highest_Number = int(str(feature.location).split("]")[0].split(":")[1])
            if int(str(feature.location).split(":")[0].strip("[")) < Lowest_Number:
                Lowest_Number = int(str(feature.location).split(":")[0].strip("["))
            #print(protein)
            if protein[2] != "0" and protein[1] != "0" and \
                    int(TM_Threshold_upper.value) >= int(protein[2]) >= int(TM_Threshold.value) and \
                    int(Cxxh_Threshold_upper.value) >= int(protein[1]) >= int(Cxxh_Threshold.value):
                color = colors.purple
            elif protein[1] != "0" and int(Cxxh_Threshold_upper.value) >= int(protein[1]) >= int(Cxxh_Threshold.value):
                color = colors.red
            elif protein[2] != "0" and int(TM_Threshold_upper.value) >= int(protein[2]) >= int(TM_Threshold.value):
                color = colors.blue

            elif len(gd_feature_set) % 2 == 0:
                color = colors.grey
            else:
                color = colors.darkgrey

            gd_feature_set.add_feature(feature, color=color, label=True)

        gd_diagram.draw(format="linear", orientation="landscape", pagesize=(int(Filesize.value)*cm, int(Filesize.value)*cm),
                        fragments=1, start=Lowest_Number, end=Highest_Number)

        gd_diagram.write("%sLinear_%s.%s" % (image_dir, file, filetype), filetype.replace(".", "")) #PS, EPS, PDF, SVG, JPG, BMP, GIF, PNG, TIFF, TIF
        # gd_diagram.write(output_dir + "Linear_" + file + "." + filetype, filetype.replace(".", "")) #PS, EPS, PDF, SVG, JPG, BMP, GIF, PNG, TIFF, TIF

        gd_diagram.draw(format="circular", circular=True, pagesize=(int(Filesize.value)*cm, int(Filesize.value)*cm),
                        start=Lowest_Number, end=Highest_Number, circle_core=0.8)
        gd_diagram.write("%sCircular_%s.%s" % (image_dir, file, filetype),
            filetype.replace(".", ""))  # PS, EPS, PDF, SVG, JPG, BMP, GIF, PNG, TIFF, TIF
    except Exception as e:
        print("Failed making figure")
        print(e)
    # gd_diagram.write(output_dir + "Circular_" + file + "." + filetype, filetype.replace(".",""))


def Spoofed_genbank(Operon, type_op, image_dir):
    try:
        global output_dir
        try:
            os.mkdir(image_dir)
        except:
            pass
        RandomInt = randrange(1000)
        Operon = Operon.split(",")
        Species = Operon[0]
        Taxid = Operon[1]
        Query = Operon[2]
        Genome = Operon[3]
        spoofed_file = open("%s%s.gb" % (image_dir, str(Query) + "_Taxid_" + Taxid + "_Genome_" + Genome + "_" + type_op), "w")  # remove random int bit ?
        #spoofed_file = open("%sFigures_TM%s_Haem%s\\%s.gb" % (output_dir, TM_Threshold.value, Cxxh_Threshold.value, str(Query) + "_Taxid_" + Taxid + "_Genome_" + Genome + "_" + type_op), "w")  # remove random int bit ?
        #spoofed_file = open(output_dir + str(Query) + "_Taxid_" + Taxid + "_Genome_" + Genome + "_" + type_op + ".gb", "w")  # remove random int bit ?
        #spoofed_file = open(output_dir + str(Query) + "_Taxid_" + Taxid +  "_" + str(RandomInt) + "_" + type_op + ".gb", "w")
        #necessary prints below
        print("LOCUS       " + Species.replace(" ", "_") + "               " + str(int(Operon[-1].split("_")[1]) - int(Operon[5].split("_")[0])) + " bp    DNA     circular BAC 53-POT-2090"
               "\n\nFEATURES             Location/Qualifiers\n     source          " + Operon[5].split("_")[0] + ".." + Operon[-1].split("_")[1], file=spoofed_file)

        #regions = []
        for i in range(4, len(Operon), 2):
            data_2 = Operon[i + 1].split("_")
            if int(data_2[-1]) < int(data_2[0]):
                data_2[-1] = str(int(data_2[0]) + int(data_2[-1]))
                Operon[i + 1] = "_".join(data_2)
            print("     gene            " + Operon[i + 1].replace("_", "..") + "\n                     /locus_tag=\"" + Operon[i] + "\"", file = spoofed_file)
            #regions.append(Operon[i + 1].split("_"))

        print("ORIGIN      ", file=spoofed_file)
        low_number = 100000000
        high_number = 0
        for checks in range(5, len(Operon), 2):
            lowest = int(Operon[checks].split("_")[0])
            highest = int(Operon[checks].split("_")[1])
            if lowest < low_number:
                low_number = lowest
            if highest > high_number:
                high_number = highest
        number = low_number  # for counting total nucleotides.
        #for i in range(1, int(Operon[-1].split("_")[1]), 60):  # start at one.
        for i in range(int(low_number), int(high_number), 60):
            first_print = "         " + str(i)
            second_print = ""
            for j in range(0, 6):
                if (number + 10) > int(Operon[-1].split("_")[1]):
                    a_left = int(Operon[-1].split("_")[1]) - number
                    second_print = second_print + "a"*a_left + "\n//"
                    #print("end")
                    break
                else:
                    second_print = second_print + "a" * 10 + " "
                    number += 10
            print(first_print[-9:] + " " + second_print, file = spoofed_file)
        spoofed_file.close()
        app.after(50, update_status, args=["Building figure...", "Yellow"])
        app.after(100, fig_maker, args=[str(Query) + "_Taxid_" + Taxid + "_Genome_" + Genome + "_" + type_op, Filetype.value, image_dir])
        #fig_maker(str(Query) + "_Taxid_" + Taxid + "_Genome_" + Genome + "_" + type_op, Filetype.value, image_dir)
        #fig_maker(str(Query) + "_" + str(RandomInt), Filetype.value)
    except Exception as e:
        print("Failed to make genbank file")
        print(e)

def add_customs():
    global customs
    customs.append(customs_value.value)
    customs_value.clear()


def windowed():
    if "custom" in ",".join(OperonType.value) or "Any" in ",".join(OperonType.value):
        window.show(wait=True)


def windowed2():
    window2.show(wait=True)


def Classifyer():
    app.after(200, update_status, args=["Building operon predictions...", "Yellow"])
    try:
        app.after(600, main_loop)
    except:
        app.after(200, update_status, args=["Failed to launch.", "Red"])

def stress_relief_update():
    #print("OK")
    if status.value.replace("-\n", "").replace(".", "") == "Building operon predictions":
        app.after(200, update_status, args=["Building operon predictions.", "Yellow"])
        app.after(200, update_status, args=["Building operon predictions..", "Yellow"])
        app.after(200, update_status, args=["Building operon predictions...", "Yellow"])
        app.after(200, update_status, args=["Building operon predictions..", "Yellow"])
        app.after(200, stress_relief_update)


#def Classifyer(TM_file, Cxxh_file, Operon_file, output_Operons, ProximateProteins, DualProteins, CxxH_threshold, TM_thresolh):
def main_loop():
    global Sequences, Haem_Weight_Dict, operon_type_dict, print_dict, customs, count_dict, top_hits, TM_thresolh, TM_thresolh_upper, \
        CxxH_threshold, CxxH_threshold_upper, output_file_name, input_csv, TM_File, cxxh_file, Haem_Weight_Ratio_To_Haem_Dict, scatterplot_data

    # print(input_csv, TM_File, cxxh_file)
    #Classifyer(TM_file, cxxh_file, input_csv, Output_File, Proximate_Only.value, Dual_Proteins.value, Cxxh_Threshold.value, TM_Threshold.value)
    #print(TM_File)
    #print(customs)
    try:
        TM_file = open(TM_File, "r")
    except:
        app.after(200, update_status, args=["ERROR- NO TM", "Red"])
        warn("No TM file!", "Please select a TM file when you import files.\nMust have \'TM\' in file name.")
        return
    try:
        Cxxh_file = open(cxxh_file, "r")
    except:
        app.after(200, update_status, args=["ERROR- NO HAEM", "Red"])
        warn("No Cxxh(Haem) file!", "Please select a CX2H/CX3H(Haem) file when you import files.\nMust have \'CX\' in file name.")
        return
    try:
        print(input_csv)
        if input_csv == "":
            p = 2 + "p"
        pass
    except:
        app.after(200, update_status, args=["ERROR- NO PROTEIN FILE", "Red"])
        warn("No protein data fie!", "Please select a protein data file when you import files.\nMust have \'Operons\' in file name.")
        return

    Sequences = False
    if "Sequences" in input_csv:
        Sequences = True
    try:
        os.mkdir("%sTM%s_%s_Haem%s_%s\\" % (str(output_dir), TM_Threshold.value,
                                            TM_Threshold_upper.value, Cxxh_Threshold.value, Cxxh_Threshold_upper.value))
    except:
        pass
    if Output_File.value == "Default":
        output_file_name = "%sTM%s_%s_Haem%s_%s\\OutputOperons%s.csv" % (str(output_dir), TM_Threshold.value, TM_Threshold_upper.value, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_"))
        output_Operons = open(output_file_name, "w")
    else:
        output_file_name = "%s%s.csv" % (str(output_dir), str(Output_File.value))
        output_Operons = open(output_file_name, "w")
    output_Operons.close()
    Operon_file = open(input_csv, "r")

    TM_dict = dict()
    Cxxh_dict = dict()
    try:
        TM_dict = {x.split("\t")[0]: x.split("\t")[1].strip() for x in TM_file}
    except:
        warn("Incorrect TM file!", "Please select a valid TM-formatted file when you import files.\nCorrect format is accession-tab-number_of_TMs")
    try:
        Cxxh_dict = {x.split("\t")[0]: x.split("\t")[1].strip() for x in Cxxh_file}
    except:
        warn("Incorrect cxxh(Haem) file!",
             "Please select a valid cxxh-formatted file when you import files.\nCorrect format is accession-tab-number_of_CXXHs")
    CxxH_threshold = int(Cxxh_Threshold.value)
    CxxH_threshold_upper = int(Cxxh_Threshold_upper.value)
    TM_thresolh = int(TM_Threshold.value)
    TM_thresolh_upper = int(TM_Threshold_upper.value)
    print_dict = dict()
    count_dict = dict()
    Haem_Weight_Dict = dict()
    Haem_Weight_Ratio_To_Haem_Dict = dict()
    scatterplot_data = dict()
    triple_count_dict = dict()
    operon_type_dict = dict()
    possible_operons = {"MtrAB": ["E-H-P-E", "E-P-H-E", "E-H-P-O", "E-P-H-O", "E-H-P-P", "E-P-H-P", "E-H-P-F",
                                  "E-P-H-F", "O-H-P-E", "O-P-H-E", "P-H-P-E", "P-P-H-E", "F-H-P-E", "F-P-H-E",
                                  "O-P-H-O", "O-P-H-P", "O-P-H-F", "P-P-H-O", "P-P-H-P", "P-P-H-F", "F-P-H-O",
                                  "F-P-H-P", "F-P-H-F", "O-H-P-O", "O-H-P-P", "O-H-P-F", "P-H-P-O", "P-H-P-P",
                                  "P-H-P-F", "F-H-P-O", "F-H-P-P", "F-H-P-F"], "MtrCAB": ["P-H-H", "H-P-H", "H-H-P"],
                        "custom": customs, "Cyc2": ["F"]}
    scatterplot_data["Other_operonname"] = []
    scatterplot_data["Other_raws"] = []
    for op_types in possible_operons:
        scatterplot_data[f"{op_types}_operonname"] = []
        scatterplot_data[f"{op_types}_raws"] = []
    try:
        oppy_namey = "_".join(OperonType.value)
        #print('%sShortData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory,OperonType.value, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value))
        saved_dict = open('%sShortData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")))
        #print("ok")
        haems_dict = open('%sHaemData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")))
        #print("ok1")
        haems_ratio_dict = open('%sHaemRatioData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")))
        #print("ok2")
        scatterplot_dict = open('%sScatterplot-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")))
        #print("ok3")
        p_dict = open('%sPrintData-%s-H%s-%s-TM%s-%s-%s.json' %
                      (sorted_variables_directory,oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")))
        #print("ok4")
        operon_type_dict = json.load(saved_dict)
        Haem_Weight_Dict = json.load(haems_dict)
        scatterplot_data = json.load(scatterplot_dict)
        Haem_Weight_Ratio_To_Haem_Dict = json.load(haems_ratio_dict)
        print_dict = json.load(p_dict)
        loaded_from_file = True
        del saved_dict
        del haems_dict
        del scatterplot_dict
        del p_dict
        #print("Loaded.")
        print(f"{datetime.now()}: All data loaded.")
    except:
        print(f"{datetime.now()}: Some data missing and needs to be generated.")
        loaded_from_file = False
        saved_dict = dict()
        haems_dict = dict()
        scatterplot_dict = dict()
        p_dict = dict()
    '''Mol_weights = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
           'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
           'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
           'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06, "_":0 }'''
    Mol_weights = {'A': 0.07104, 'C': 0.10301, 'D': 0.11503, 'E': 0.12904, 'F': 0.14707,
           'G': 0.05702, 'H': 0.13706, 'I': 0.11308, 'K': 0.12809, 'L': 0.11308,
           'M': 0.13104, 'N': 0.11404, 'P': 0.09705, 'Q': 0.12806, 'R': 0.15610,
           'S': 0.08703, 'T': 0.10105, 'V': 0.09907, 'W': 0.18608, 'Y': 0.16306, "_": 0}

    #weight = sum(Mol_weights[p] for p in seq)

    print("Operon Type,TaxID,Query Sequence,Whole Genome Sequenced Match,ProteinID-Cxxh_motifs-TM_regions-MolWeight----->", file = open(output_file_name, "w"))
    # print("Settings\nSequences:%s\nLoaded from files:%s" % (Sequences, loaded_from_file))
    if Sequences and not loaded_from_file:
        for line in Operon_file:
            data = line.strip().split(",")
            #query = data[1]  # old one... could have many hits! this only displays one!
            query = []
            #print(data)
            if hide_query_duplicate.value == "Hide duplicate query ID":  # make this togglable option
                query = [data[0], data[2]]  # displays unique to taxid+target genome + hit
                for prot in range(3, len(data), 2):
                    query.append(data[prot].replace("*", ""))
            else:
                query = [data[0], data[1], data[2]]  # taxid + query + genome + operon
                for prot in range(3, len(data), 2):
                    query.append(data[prot].replace("*", ""))
                # query = ",".join([data[0], data[1], data[2], data[3]])  # displays unique to taxid+query + target genome
            Out_data = []

            query = ",".join(query)
            #print(query)
            count_dict[query] = 0
            operon_descriptor = ["E"]
            Haem_Weight_Dict[query] = 0
            Haem_Weight_Ratio_To_Haem_Dict[query] = []  # list as unique numbers can repat. overlap won't show on scatterplot, but i dont like to lose data
            operon_type_dict[query] = []
            print_dict[query] = ",".join([data[0], data[1], data[2]])
            for protein in range(3, len(data), 2):
                #cell = data[protein].strip().split("(")[0].strip("*")
                accession = data[protein].strip().split("(")[0].strip("*")
                cell1 = data[protein].strip()
                sequence = data[protein + 1].strip().strip("*")
                weight = sum(Mol_weights[p] for p in sequence)

                if weight == 0:
                    weight = 10000000

                if accession in TM_dict and accession in Cxxh_dict:
                    print_dict[query] = print_dict[query] + "," + cell1 + "-" + \
                                        str(int(Cxxh_dict[accession.strip("*")])) + "-" + str(int(TM_dict[accession.strip("*")])) + "-" \
                                        + str(weight)
                    if int(TM_dict[accession.strip("*")]) >= TM_thresolh and int(TM_dict[accession.strip("*")]) <= TM_thresolh_upper and int(Cxxh_dict[accession.strip("*")]) >= CxxH_threshold and int(Cxxh_dict[accession.strip("*")]) <= CxxH_threshold_upper:
                        operon_descriptor.append("F")
                        haem_number = float(Cxxh_dict[accession.strip("*")])
                        if float(Haem_Weight_Dict[query]) < float(weight/haem_number):
                            Haem_Weight_Dict[query] = round(float(weight/haem_number), 1)
                        Haem_Weight_Ratio_To_Haem_Dict[query].append(f"{int(haem_number)}:{round(float(weight/haem_number), 1)}")

                    elif int(Cxxh_dict[accession.strip("*")]) >= CxxH_threshold and int(Cxxh_dict[accession.strip("*")]) <= CxxH_threshold_upper:
                        operon_descriptor.append("H")
                        haem_number = float(Cxxh_dict[accession.strip("*")])
                        if float(Haem_Weight_Dict[query]) < float(weight/haem_number):
                            Haem_Weight_Dict[query] = round(float(weight/haem_number), 1)
                        Haem_Weight_Ratio_To_Haem_Dict[query].append(f"{int(haem_number)}:{round(float(weight/haem_number), 1)}")
                    elif int(TM_dict[accession.strip("*")]) >= TM_thresolh and int(TM_dict[accession.strip("*")]) <= TM_thresolh_upper:
                        operon_descriptor.append("P")
                    elif int(Cxxh_dict[accession.strip("*")]) > 0:
                        operon_descriptor.append("h")
                    else:
                        operon_descriptor.append("O")


                elif accession in TM_dict:
                    print_dict[query] = print_dict[query] + "," + cell1 + "-0-" + str(int(TM_dict[accession.strip("*")]) ) + "-" + str(weight)
                    if int(TM_dict[accession.strip("*")]) >= TM_thresolh and int(TM_dict[accession.strip("*")]) <= TM_thresolh_upper:
                        operon_descriptor.append("P")
                    else:
                        operon_descriptor.append("O")

                elif accession in Cxxh_dict:
                    print_dict[query] = print_dict[query] + "," + cell1 + "-" + \
                                        str(int(Cxxh_dict[accession.strip("*")])) + "-0-" + str(weight)
                    if int(Cxxh_dict[accession.strip("*")]) >= CxxH_threshold and int(Cxxh_dict[accession.strip("*")]) <= CxxH_threshold_upper:
                        operon_descriptor.append("H")
                        haem_number = float(Cxxh_dict[accession.strip("*")])
                        if float(Haem_Weight_Dict[query]) < float(weight/haem_number):
                            Haem_Weight_Dict[query] = round(float(weight/haem_number), 1)
                        Haem_Weight_Ratio_To_Haem_Dict[query].append(f"{int(haem_number)}:{round(float(weight/haem_number), 1)}")
                    elif int(Cxxh_dict[accession.strip("*")]) > 0:
                        operon_descriptor.append("h")
                        haem_number = float(Cxxh_dict[accession.strip("*")])
                        if float(Haem_Weight_Dict[query]) < float(weight/haem_number):
                            Haem_Weight_Dict[query] = round(float(weight/haem_number), 1)
                        Haem_Weight_Ratio_To_Haem_Dict[query].append(f"{int(haem_number)}:{round(float(weight/haem_number), 1)}")
                    else:
                        operon_descriptor.append("O")
                else:
                    print_dict[query] = print_dict[query] + "," + cell1 + "-0-0-" + str(weight)
                    operon_descriptor.append("O")
            ##
            operon_descriptor.append("E")
            operon_descriptor = "-".join(operon_descriptor)

            #print(operon_descriptor)
            if "Any" in ".".join(OperonType.value):
                for key in possible_operons:
                    for op_type in possible_operons[key]:
                        if op_type in operon_descriptor.replace("h", "O"):
                            operon_type_dict[query].append(key)
                            scatterplot_data[f"{key}_operonname"].append(query)
                            scatterplot_data[f"{key}_raws"].append(Haem_Weight_Ratio_To_Haem_Dict[query])
                            #scatterplot_data[key] = scatterplot_data[key] + Haem_Weight_Ratio_To_Haem_Dict[query]
                if not operon_type_dict[query]:
                    if ("P" in operon_descriptor.replace("F", "P") and "H" in operon_descriptor) or \
                            ("H" in operon_descriptor.replace("F", "H") and "P" in operon_descriptor):
                        operon_type_dict[query].append("Other")
                        scatterplot_data["Other_operonname"].append(query)
                        scatterplot_data["Other_raws"].append(Haem_Weight_Ratio_To_Haem_Dict[query])
            else:
                for operon_type in OperonType.value:
                    for op_type in possible_operons[operon_type]:
                        if op_type in operon_descriptor:
                            operon_type_dict[query].append(operon_type)
                            scatterplot_data[operon_type + "_operonname"].append(query)
                            scatterplot_data[operon_type + "_raws"].append(Haem_Weight_Ratio_To_Haem_Dict[query])
            """if OperonType.value == "Cyc2":
                for op_type in possible_operons["Cyc2"]:
                    if op_type in operon_descriptor:
                        operon_type_dict[query].append("Cyc2")
                        scatterplot_data["Cyc2_operonname"].append(query)
                        scatterplot_data["Cyc2_raws"].append(Haem_Weight_Ratio_To_Haem_Dict[query])
            elif OperonType.value == "MtrAB":
                for op_type in possible_operons["MtrAB"]:
                    if op_type in operon_descriptor:
                        operon_type_dict[query].append("MtrAB")
                        scatterplot_data["MtrAB_operonname"].append(query)
                        scatterplot_data["MtrAB_raws"].append(Haem_Weight_Ratio_To_Haem_Dict[query])
            elif OperonType.value == "MtrCAB":
                for op_type in possible_operons["MtrCAB"]:
                    if op_type in operon_descriptor:
                        operon_type_dict[query].append("MtrCAB")
                        scatterplot_data["MtrCAB_operonname"].append(query)
                        scatterplot_data["MtrCAB_raws"].append(Haem_Weight_Ratio_To_Haem_Dict[query])
            elif OperonType.value == "custom":
                for op_type in possible_operons["custom"]:
                    if op_type in operon_descriptor:
                        operon_type_dict[query].append("custom")
                        scatterplot_data["custom_operonname"].append(query)
                        scatterplot_data["custom_raws"].append(Haem_Weight_Ratio_To_Haem_Dict[query])"""
            operon_type_dict[query] = ";".join(set(operon_type_dict[query]))

            #print(operon_type_dict[query])
        ### ENDS HERE
    elif not loaded_from_file:
        for line in Operon_file:
            #data = line.strip().split("(")[0].split(",")
            data = line.strip().split(",")
            if hide_query_duplicate.value == "Hide duplicate query ID":  # make this togglable option
                query = [data[0]]  # displays unique to taxid+target genome + hit
                for prot in range(3, len(data), 2):
                    query.append(data[prot].replace("*", ""))
            else:
                query = [data[0], data[1], data[2]]  # taxid + query + genome + operon
                for prot in range(3, len(data)): #no sequences, so no 2-step
                    query.append(data[prot].replace("*", ""))
            Out_data = []
            query = ",".join(query)
            count_dict[query] = 0
            operon_descriptor = ["E"]
            operon_type_dict[query] = ""
            print_dict[query] = ",".join([data[0], data[1], data[2]])
            for protein in range(3, len(data)):

                accession = data[protein].strip().split("(")[0].strip("*")
                cell1 = data[protein]
                if accession in TM_dict and accession in Cxxh_dict:
                    print_dict[query] = print_dict[query] + "," + cell1 + "-" + \
                                        str(int(Cxxh_dict[accession.strip("*")])) + "-" + str(int(TM_dict[accession.strip("*")]))
                    if int(TM_dict[accession.strip("*")]) >= TM_thresolh and int(Cxxh_dict[accession.strip("*")]) >= CxxH_threshold\
                            and int(TM_dict[accession.strip("*")]) <= TM_thresolh_upper and int(Cxxh_dict[accession.strip("*")]) <= CxxH_threshold_upper:
                        operon_descriptor.append("F")
                        if query in count_dict:
                            count_dict[query] = count_dict[query] + 2
                        else:
                            count_dict[query] = 2
                    elif int(Cxxh_dict[accession.strip("*")]) >= CxxH_threshold and int(Cxxh_dict[accession.strip("*")]) <= CxxH_threshold_upper:
                        operon_descriptor.append("H")
                        if query in count_dict:
                            count_dict[query] = count_dict[query] + 1
                        else:
                            count_dict[query] = 1
                    elif int(TM_dict[accession.strip("*")]) >= TM_thresolh and int(TM_dict[accession.strip("*")]) <= TM_thresolh_upper:
                        operon_descriptor.append("P")
                        if query in count_dict:
                            count_dict[query] = count_dict[query] + 1
                        else:
                            count_dict[query] = 1
                    elif int(Cxxh_dict[accession.strip("*")]) >= 0:
                        operon_descriptor.append("h")
                    else:
                        operon_descriptor.append("O")

                elif accession in TM_dict:
                    print_dict[query] = print_dict[query] + "," + cell1 + "-0-" + str(int(TM_dict[accession.strip("*")]))
                    if int(TM_dict[accession.strip("*")]) >= TM_thresolh and int(TM_dict[accession.strip("*")]) <= TM_thresolh_upper:
                        operon_descriptor.append("P")
                        if query in count_dict:
                            count_dict[query] = count_dict[query] + 1
                        else:
                            count_dict[query] = 1
                    else:
                        operon_descriptor.append("O")
                elif accession in Cxxh_dict:
                    print_dict[query] = print_dict[query] + "," + cell1 + "-" + \
                                        str(int(Cxxh_dict[accession.strip("*")])) + "-0"
                    if int(Cxxh_dict[accession.strip("*")]) >= CxxH_threshold and int(Cxxh_dict[accession.strip("*")]) <= CxxH_threshold_upper:
                        operon_descriptor.append("H")
                        if query in count_dict:
                            count_dict[query] = count_dict[query] + 1
                        else:
                            count_dict[query] = 1
                    elif int(Cxxh_dict[accession.strip("*")]) >= 0:
                        operon_descriptor.append("h")
                    else:
                        operon_descriptor.append("O")
                else:
                    print_dict[query] = print_dict[query] + "," + cell1 + "-0-0"
                    operon_descriptor.append("O")
            operon_descriptor.append("E")
            operon_descriptor = "-".join(operon_descriptor)
            if "Any" in ".".join(OperonType.value):
                for key in possible_operons:
                    for op_type in possible_operons[key]:
                        if op_type in operon_descriptor.replace("h", "O"):
                            operon_type_dict[query] = operon_type_dict[query] + [key]
                if not operon_type_dict[query]:
                    if ("P" in operon_descriptor.replace("F", "P") and "H" in operon_descriptor) or \
                            ("H" in operon_descriptor.replace("F", "H") and "P" in operon_descriptor):
                        operon_type_dict[query] = operon_type_dict[query] + ["Other"]
            else:
                for operon_type in OperonType.value:
                    for op_type in possible_operons[operon_type]:
                        if op_type in operon_descriptor:
                            operon_type_dict[query] = operon_type_dict[query] + [operon_type]
            """if OperonType.value == "Cyc2":
                for op_type in possible_operons["Cyc2"]:
                    if op_type in operon_descriptor:
                        operon_type_dict[query] = operon_type_dict[query] + ["Cyc2"]
            elif OperonType.value == "MtrAB":
                for op_type in possible_operons["MtrAB"]:
                    if op_type in operon_descriptor:
                        operon_type_dict[query] = operon_type_dict[query] + ["MtrAB"]
            elif OperonType.value == "MtrCAB":
                for op_type in possible_operons["MtrCAB"]:
                    if op_type in operon_descriptor:
                        operon_type_dict[query] = operon_type_dict[query] + ["MtrCAB"]
            elif OperonType.value == "custom":
                for op_type in possible_operons["custom"]:
                    if op_type in operon_descriptor:
                        operon_type_dict[query] = operon_type_dict[query] + ["custom"]"""
            operon_type_dict[query] = ";".join(set(operon_type_dict[query]))
    top_hits = int(Top_Hits.value)
    oppy_namey = "_".join(OperonType.value)
    if not os.path.exists('%sShortData-%s-H%s-%s-TM%s-%s-%s.json' %
                          (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_"))):
        with open('%sShortData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")), 'w') as f:
            saved_dict = operon_type_dict
            json.dump(saved_dict, f)
            del saved_dict
    if not os.path.exists('%sHaemData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_"))):
        with open('%sHaemData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")), 'w') as f:
            haems_dict = Haem_Weight_Dict
            json.dump(haems_dict, f)
            del haems_dict
    if not os.path.exists('%sScatterplot-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_"))):
        with open('%sScatterplot-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")), 'w') as f:
            scatterplot_dict = scatterplot_data
            json.dump(scatterplot_dict, f)
            del scatterplot_dict

    if not os.path.exists('%sHaemRatioData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_"))):
        with open('%sHaemRatioData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")), 'w') as f:
            haems_ratio_dict = Haem_Weight_Ratio_To_Haem_Dict
            json.dump(haems_ratio_dict, f)
            del haems_ratio_dict

    if not os.path.exists('%sPrintData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_"))):
        with open('%sPrintData-%s-H%s-%s-TM%s-%s-%s.json' % (sorted_variables_directory, oppy_namey, Cxxh_Threshold.value, Cxxh_Threshold_upper.value, TM_Threshold.value, TM_Threshold_upper.value, hide_query_duplicate.value.replace(" ", "_")), 'w') as f:
            p_dict = print_dict
            json.dump(p_dict, f)
            del p_dict
    print(f"{datetime.now()}: Data all saved and ready to create output figures.")
    app.after(500, update_status, args=["Data analysed...", "Orange"])
    app.after(1000, update_status, args=["Outputting data...", "Orange"])
    try:
        app.after(1500, last_bit)
    except:
        app.after(1500, update_status, args=["ERROR- Failed data output", "Red"])


def update_phy_and_itol(phy_tree):
    try:
        if phy_tree.endswith(".phy"):
            all_lines = []
            for line in open(phy_tree, "r"):
                all_lines.append(line.replace("\':", "==tempy==").replace(":", "_").replace("==tempy==", "\':")
                                 .replace(",\n", "==tempy==").replace(",", "_").replace("==tempy==", ",\n").strip())
            print("\n".join(all_lines), file=open(phy_tree.replace(".phy", "_Updated.phy"), "w"))
            app.after(800, update_status, args=[f"Updated iTOL phy tree!", "Green"])
        else:
            app.after(800, update_status, args=[f"ERROR Updating phy tree (Not .phy file)!", "Red"])
    except Exception as e:
        app.after(800, update_status, args=[f"ERROR Updating phy tree ({e})!", "Red"])


def make_itol_files(output_dir, TM_thresolh, CxxH_threshold, histogram_data, TM_thresolh_upper, CxxH_threshold_upper, scatterplot_data, operon_type_dict):
    # operon_type_dict = {QUERY: ["Cyc2"], QUERY2: ["MtrAB", "Cyc2"]...}
    #print(operon_type_dict)
    print(f"{datetime.now()}: Started on iTOL files.")
    #time.sleep(30)
    itol_output = "%sTM%s_%s_Haem%s_%s\\iTOL_Files\\" % \
                  (output_dir, TM_Threshold.value, TM_Threshold_upper.value, Cxxh_Threshold.value, Cxxh_Threshold_upper.value)
    try:
        os.mkdir(itol_output)
    except:
        pass
    if itol_many_genome_files.value == "All genome files (iTOL)":
        generic_filename = "%sGeneric_iTOL.txt" % (itol_output)
    else:
        generic_filename = "%sGeneric_iTOL_TopGenomeOnly.txt" % (itol_output)
    generic_itol_file = open(generic_filename, "w")
    generic_itol_file.close()
    if itol_many_genome_files.value == "All genome files (iTOL)":
        specific_filename = "%sSpecific_iTOL.txt" % (itol_output)
    else:
        specific_filename = "%sSpecific_iTOL_TopGenomeOnly.txt" % (itol_output)
    specific_itol_file = open(specific_filename, "w")
    specific_itol_file.close()

    instructions = ["Using file TaxonList.txt", "",
                    "Go to : https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi", "",
                    "Upload the file of taxonomic IDs and click add from file", "",
                    "Click the check box next to root to select everything and click choose", "",
                    "Click the drop down box next to save as and select phylip tree, then click save as", "",
                    "Make sure to save the file here! ", "", "Run the program again with the same settings.",
                    "If you get the \'Failed Name Conversion\' warning, you will have to re-download names.dmp "
                    "file from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/.",
                    "Download https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip and unzip it into the home directory",
                    "Replace the taxdmp folder there.",
                    "", "You will probably want the option \'One Genome File(iTOL)\' that ensures only one genome file "
                        "is used per Taxid (one with most predicted operons). ",
                    "\'Hide Duplicate Query ID\' is also very useful, as you remove duplicates from multiple query "
                    "sequences matching a single hit in the database.", "", "Making tree:",
                    "1) Load the phylip tree file into itol to load the tree.", "2) Add Specific_iTOL.txt or"
                                                                                " Generic_iTOL.txt to dataset"]
    print("\n".join(instructions), file=open("%sInstructions.txt" % (itol_output), "w"))
    clear_1 = open("%sTaxon_to_Species.txt" % (itol_output), "w")
    clear_2 = open("%sMissing_Taxons.txt" % (itol_output), "w")
    clear_3 = open("%sTaxonList.txt" % (itol_output), "w")
    clear_1.close()
    clear_2.close()
    clear_3.close()
    changed_taxids = {}
    backup_taxids = {}
    for line in open(currdir + "taxdmp\\merged.dmp", "r"):
        data = line.split("|")
        changed_taxids[data[0].strip().replace("\t", "")] = data[1].strip().replace("\t", "")
        backup_taxids[data[1].strip().replace("\t", "")] = data[0].strip().replace("\t", "")

    taxid_to_name_dict = {}
    speciesname_to_taxid_dict = {}
    blank_itol_output_dict = {}
    for line in open(currdir + "taxdmp\\names.dmp", "r"):
        # newest from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
        # get new file if there are errors
        data = line.split("|")  # taxid\t|\tspecies name
        taxid_here = data[0].strip().replace("\t", "")
        species_here = data[1].strip().replace("\t", "").replace(":", "_").replace(",", "_").replace("\'", "")
        if taxid_here in changed_taxids: # if there is an old and new version, change just in case
            taxid_to_name_dict[changed_taxids[taxid_here]] = species_here
            taxid_to_name_dict[taxid_here] = species_here
            speciesname_to_taxid_dict[species_here] = \
                changed_taxids[taxid_here]
            if itol_many_genome_files.value == "All genome files (iTOL)":
                blank_itol_output_dict[taxid_here] = {"Cyc2": 0, "MtrAB": 0, "MtrCAB": 0, "Other": 0, "custom": 0,
                                                "SpeciesName": species_here, "Total": 0}
            else:
                blank_itol_output_dict[taxid_here] = {
                    "Blank": {"Cyc2": 0, "MtrAB": 0, "MtrCAB": 0, "Other": 0, "custom": 0,
                              "SpeciesName": species_here, "Total": 0}}

        else:
            taxid_to_name_dict[taxid_here] = species_here
            speciesname_to_taxid_dict[species_here] = taxid_here
            if itol_many_genome_files.value == "All genome files (iTOL)":
                blank_itol_output_dict[taxid_here] = {"Cyc2": 0, "MtrAB": 0, "MtrCAB": 0, "Other": 0, "custom": 0,
                                                   "SpeciesName": species_here, "Total": 0}
            else:
                blank_itol_output_dict[taxid_here] = {
                    "Blank": {"Cyc2": 0, "MtrAB": 0, "MtrCAB": 0, "Other": 0, "custom": 0,
                              "SpeciesName": species_here, "Total": 0}}


    itol_output_dict = {}
     # itol_many_genome_files.value == All genome files (iTOL)
    operon_types_present_list = []
    if itol_many_genome_files.value == "All genome files (iTOL)":
        for each_query in operon_type_dict:
            operon_taxid = each_query.split(",")[0].strip(")").split("(")[-1]  # tax IDs like: 651,383219,398137,..
            if hide_query_duplicate.value == "Do not hide":
                current_genome = each_query.split(",")[2]
            else:
                current_genome = each_query.split(",")[1]
            if operon_taxid not in blank_itol_output_dict:
                if operon_taxid in changed_taxids:
                    sp_name = taxid_to_name_dict[changed_taxids[operon_taxid]]  # if converting old taxid, do it
                elif operon_taxid in taxid_to_name_dict:
                    sp_name = taxid_to_name_dict[operon_taxid]  # if correct taxid, change to species name
                else:
                    sp_name = operon_taxid  # if no species name, use taxid for now...
                blank_itol_output_dict[operon_taxid] = {"Cyc2": 0, "MtrAB": 0, "MtrCAB": 0, "Other": 0, "custom": 0,
                                                        "SpeciesName": sp_name.replace(":", "_").replace(",", "_"),
                                                        "Total": 0}

            if operon_taxid in blank_itol_output_dict and operon_type_dict[each_query] != "":
                for operon_type in operon_type_dict[each_query].split(";"):
                    operon_types_present_list.append(operon_type)
                    blank_itol_output_dict[operon_taxid][operon_type] += 1
                    blank_itol_output_dict[operon_taxid]["Total"] += 1
    else:
        for each_query in operon_type_dict:
            operon_taxid = each_query.split(",")[0].strip(")").split("(")[-1]  # tax IDs like: 651,383219,398137,..
            if hide_query_duplicate.value == "Do not hide":
                current_genome = each_query.split(",")[2]
            else:
                current_genome = each_query.split(",")[1]
            if operon_taxid not in blank_itol_output_dict:
                if operon_taxid in changed_taxids:
                    sp_name = taxid_to_name_dict[changed_taxids[operon_taxid]]  # if converting old taxid, do it
                elif operon_taxid in taxid_to_name_dict:
                    sp_name = taxid_to_name_dict[operon_taxid]  # if correct taxid, change to species name
                else:
                    sp_name = operon_taxid  # if no species name, use taxid for now...
                blank_itol_output_dict[operon_taxid] = {current_genome: {"Cyc2": 0, "MtrAB": 0, "MtrCAB": 0, "Other": 0, "custom": 0,
                                                        "SpeciesName": sp_name.replace(":", "_").replace(",", "_"), "Total": 0}}
            elif current_genome not in blank_itol_output_dict[operon_taxid]:
                if operon_taxid in changed_taxids:
                    sp_name = taxid_to_name_dict[changed_taxids[operon_taxid]]  # if converting old taxid, do it
                elif operon_taxid in taxid_to_name_dict:
                    sp_name = taxid_to_name_dict[operon_taxid]  # if correct taxid, change to species name
                else:
                    sp_name = operon_taxid  # if no species name, use taxid for now...
                blank_itol_output_dict[operon_taxid][current_genome] = {"Cyc2": 0, "MtrAB": 0, "MtrCAB": 0,
                                                                        "Other": 0, "custom": 0, "Total": 0,
                                                                        "SpeciesName": sp_name.replace(":", "_").replace(",", "_")}
            if operon_taxid in blank_itol_output_dict and operon_type_dict[each_query] != "":
                for operon_type in operon_type_dict[each_query].split(";"):
                    operon_types_present_list.append(operon_type)
                    blank_itol_output_dict[operon_taxid][current_genome][operon_type] += 1
                    blank_itol_output_dict[operon_taxid][current_genome]["Total"] += 1
        top_genome_list = []
        for each_taxid in blank_itol_output_dict:
            top_genome = "None yet"
            top_genome_score = -1
            for each_genome in blank_itol_output_dict[each_taxid]:
                current_minidict = blank_itol_output_dict[each_taxid][each_genome]
                if (current_minidict["Total"]) > top_genome_score:
                    top_genome = each_genome
                    top_genome_score = current_minidict["Total"]
            top_genome_list.append(f"{each_taxid}:{top_genome}")
            blank_itol_output_dict[each_taxid] = blank_itol_output_dict[each_taxid][top_genome]
        print("\n".join(top_genome_list), file=open("%sTaxID_TopGenome_List.txt" % (itol_output), "w"))
    operon_types_present_list = list(set(operon_types_present_list))
    colours_for_itol = ["#f42b1d", "#0077e1", "#BADA55", "#a13fff", "#ffa03d"][0:len(operon_types_present_list)]
    itol_output_generic = ["DATASET_MULTIBAR",
                           "#In multi-value bar charts, each ID is associated to multiple numeric values, which are displayed as a stacked or aligned bar chart",
                           "#lines starting with a hash are comments and ignored during parsing", "",
                           "#=================================================================#",
                           "#                    MANDATORY SETTINGS                           #",
                           "#=================================================================#",
                           "#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.",
                           "#SEPARATOR TAB", "#SEPARATOR SPACE", "SEPARATOR COMMA", "",
                           "#label is used in the legend table (can be changed later)",
                           "DATASET_LABEL,example multi bar chart", "", "#dataset color (can be changed later)",
                           "COLOR,#ff0000", "",
                           "#define colors for each individual field column (use hexadecimal, RGB or RGBA notation; if using RGB/RGBA, COMMA cannot be used as SEPARATOR)",
                           "FIELD_COLORS,#ff0000,#00ff00,#0000ff", "", "#field labels", "FIELD_LABELS,f1,f2,f3", "",
                           "#=================================================================#",
                           "#                    OPTIONAL SETTINGS                            #",
                           "#=================================================================#", "", "", "",
                           "#=================================================================#",
                           "#     all other optional settings can be set or changed later     #",
                           "#           in the web interface (under 'Datasets' tab)           #",
                           "#=================================================================#", "",
                           "#dataset scale: you can simply set the values where the scale will be drawn",
                           "#DATASET_SCALE,2000,10000,20000",
                           "#or you can specify value, label, color, width, style and label size factor for each scale line (dash separated, format: VALUE-LABEL-COLOR-WIDTH-DASHED-LABEL_SCALE_FACTOR))",
                           "#DATASET_SCALE,2000-2k line-#0000ff-5-1-1,10000-line at 10k-#ff0000-1-0-2,20000-3rd line-#00ff00-5-1-1",
                           "", "", "#Each dataset can have a legend, which is defined using LEGEND_XXX fields below",
                           "#For each row in the legend, there should be one shape, color and label.",
                           "#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values",
                           "#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.",
                           "#Shape should be a number between 1 and 6, or any protein domain shape definition.",
                           "#1: square", "#2: circle", "#3: star", "#4: right pointing triangle",
                           "#5: left pointing triangle", "#6: checkmark", "", "#LEGEND_TITLE,Dataset legend",
                           "#LEGEND_POSITION_X,100", "#LEGEND_POSITION_Y,100", "#LEGEND_SHAPES,1,2,3",
                           "#LEGEND_COLORS,#ff0000,#00ff00,#0000ff", "#LEGEND_LABELS,value1,value2,value3",
                           "#LEGEND_SHAPE_SCALES,1,1,0.5", "", "#maximum width", "#WIDTH,1000", "",
                           "#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.",
                           "#MARGIN,0", "",
                           "#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.",
                           "#SHOW_INTERNAL,0", "", "#show dashed lines between leaf labels and the dataset",
                           "DASHED_LINES,1", "",
                           "#bar height factor; Default bar height will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it, values above 1 will increase it)",
                           "#HEIGHT_FACTOR,1", "",
                           "#Bars are aligned to the node lines by default. Using BAR_SHIFT, you can move them all up/down by a fixed amount",
                           "#BAR_SHIFT,0", "",
                           "#align individual fields; if set to 1, individual bar charts will not be stacked",
                           "#ALIGN_FIELDS,0", "",
                           "#border width; if set above 0, a border of specified width (in pixels) will be drawn around the bars",
                           "#BORDER_WIDTH,0", "", "#border color; used when BORDER_WIDTH is above 0",
                           "#BORDER_COLOR,#0000ff", "",
                           "#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages",
                           "#=================================================================#",
                           "#       Actual data follows after the \"DATA\" keyword              #",
                           "#=================================================================#", "DATA"]
    itol_output_readytogo = ["DATASET_MULTIBAR", "", "SEPARATOR COMMA", "",
                             "#label is used in the legend table (can be changed later)",
                             "DATASET_LABEL, ETMiner Multibar Chart", "", "#dataset color (can be changed later)",
                             "COLOR,#2a009d", "",
                             "#define colors for each individual field column (use hexadecimal, RGB or RGBA notation; if using RGB/RGBA, COMMA cannot be used as SEPARATOR)",
                             "FIELD_COLORS," + ",".join(colours_for_itol), "", "#field labels",
                             "FIELD_LABELS," + ",".join(operon_types_present_list), "", "",
                             "#dataset scale: you can simply set the values where the scale will be drawn",
                             "#DATASET_SCALE,2000,10000,20000",
                             "#or you can specify value, label, color, width, style and label size factor for each scale line (dash separated, format: VALUE-LABEL-COLOR-WIDTH-DASHED-LABEL_SCALE_FACTOR))",
                             "#DATASET_SCALE,2000-2k line-#0000ff-5-1-1,10000-line at 10k-#ff0000-1-0-2,20000-3rd line-#00ff00-5-1-1",
                             "", "", "#Each dataset can have a legend, which is defined using LEGEND_XXX fields below",
                             "#For each row in the legend, there should be one shape, color and label.",
                             "#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values",
                             "#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.",
                             "#Shape should be a number between 1 and 6, or any protein domain shape definition.",
                             "#1: square", "#2: circle", "#3: star", "#4: right pointing triangle",
                             "#5: left pointing triangle", "#6: checkmark", "", "LEGEND_TITLE, Legend",
                             "#LEGEND_POSITION_X,100", "#LEGEND_POSITION_Y,100", "LEGEND_SHAPES," + ",".join(["4"]*len(operon_types_present_list)),
                             "LEGEND_COLORS," + ",".join(colours_for_itol),
                             "LEGEND_LABELS," + ",".join(operon_types_present_list), "#LEGEND_SHAPE_SCALES,1,1,0.5", "", "DATA"]
    all_warns = []
    for phy_tree in glob.iglob("%s*.phy" % itol_output):
        if "updated" not in phy_tree.lower():
            all_lines = []
            for line in open(phy_tree, "r"):
                bit_to_replace = "made up placeholder"
                bit_to_replace_temp = "made up placeholder"
                newline = line
                if "\'" in line:
                    newline = line.replace("\':", "==tempy==").replace(":", "_").replace("==tempy==", "\':")\
                        .replace(",\n", "==tempy==").replace(",", "_").replace("==tempy==", ",\n").strip()
                # newline == "\'NCBI SPECIES NAME\':leaf info"
                    bit_to_replace = newline.split(":")[0]
                    bit_to_replace_temp = bit_to_replace.strip().replace("\'", "").replace(";", "").replace(",", "")
                if bit_to_replace_temp[0] == ")":
                    bit_to_replace_temp = bit_to_replace_temp[1:]
                if bit_to_replace[0] == ")":
                    bit_to_replace = bit_to_replace[1:]
                if bit_to_replace_temp in speciesname_to_taxid_dict :  # if species name is in the .phy file, set it as the main for output
                    if speciesname_to_taxid_dict[bit_to_replace_temp] in blank_itol_output_dict:
                        #print(bit_to_replace)
                        blank_itol_output_dict[speciesname_to_taxid_dict[bit_to_replace_temp]]["SpeciesName"] = bit_to_replace_temp.replace("(", "").replace(")", "").replace("[", "").replace("]", "")
                        newline = newline.replace(bit_to_replace, "\'" + blank_itol_output_dict[speciesname_to_taxid_dict[bit_to_replace_temp]]["SpeciesName"] + "\'")
                        itol_output_dict[speciesname_to_taxid_dict[bit_to_replace_temp]] = blank_itol_output_dict[speciesname_to_taxid_dict[bit_to_replace_temp]]
                elif bit_to_replace_temp != "made up placeholder":
                    all_warns.append(f"Failed to replace \'{bit_to_replace_temp}\' from line \'{newline}\'")
                    print(f"Failed to replace \'{bit_to_replace_temp}\' from line \'{newline}\' -- Check names.dmp"
                          f" file!")
                    if bit_to_replace_temp in speciesname_to_taxid_dict:
                        print(f"\'{bit_to_replace_temp}\' is in taxid system. Taxid:{speciesname_to_taxid_dict[bit_to_replace_temp]}"
                              f"\nAdded to output in mock format.")
                    itol_output_dict[bit_to_replace_temp] = {"Cyc2": 0, "MtrAB": 0, "MtrCAB": 0, "Other": 0, "custom": 0,
                                                        "SpeciesName": bit_to_replace_temp}

                all_lines.append(newline.strip())
            if all_warns:
                all_warns = "\n".join(all_warns)
                warn("Failed Name Conversion", f"Failed on following species-names-taxid conversions:\n{all_warns}")
            print("\n".join(all_lines), file=open(phy_tree.replace(".phy", "_Updated.phy"), "w"))
    #print("\n".join(list(itol_output_dict.keys())), file=open("%sTaxonList.txt" % itol_output, "w"))
    if list(itol_output_dict.keys()):
        print("\n".join(list(itol_output_dict.keys())), file=open("%sTaxonList.txt" % itol_output, "w"))
    else:
        send_list = [x for x in blank_itol_output_dict if blank_itol_output_dict[x]["Total"] > 0]
        print("\n".join(send_list), file=open("%sTaxonList.txt" % itol_output, "w"))
    taxa_for_file = []
    for species_taxon_id in itol_output_dict:  # Cyc2,MtrAB,MtrCAB,Other,Custom
        taxa_for_file.append(species_taxon_id)
        for_itol_output_mini = ["\'" + str(itol_output_dict[species_taxon_id]["SpeciesName"]) + "\'"]
        for op_typer in operon_types_present_list:
            for_itol_output_mini.append(str(itol_output_dict[species_taxon_id][op_typer]))
        itol_output_generic.append(",".join(for_itol_output_mini))
        itol_output_readytogo.append(",".join(for_itol_output_mini))
        """itol_output_generic.append(",".join(["\'" + str(itol_output_dict[species_taxon_id]["SpeciesName"]) + "\'",
                                             str(itol_output_dict[species_taxon_id]["Cyc2"]),
                                             str(itol_output_dict[species_taxon_id]["MtrAB"]),
                                             str(itol_output_dict[species_taxon_id]["MtrCAB"]),
                                             str(itol_output_dict[species_taxon_id]["Other"]),
                                             str(itol_output_dict[species_taxon_id]["custom"])]))
        itol_output_readytogo.append(",".join(["\'" + str(itol_output_dict[species_taxon_id]["SpeciesName"]) + "\'",
                                               str(itol_output_dict[species_taxon_id]["Cyc2"]),
                                               str(itol_output_dict[species_taxon_id]["MtrAB"]),
                                               str(itol_output_dict[species_taxon_id]["MtrCAB"]),
                                               str(itol_output_dict[species_taxon_id]["Other"]),
                                               str(itol_output_dict[species_taxon_id]["custom"])]))"""
    all_taxons_print = []
    all_missing_taxons_print = []
    for each_taxid in taxa_for_file:
        if each_taxid in taxid_to_name_dict:
            all_taxons_print.append(":".join([each_taxid, taxid_to_name_dict[each_taxid].replace(":", "_").replace(",", "_").replace("\'", "")]))
        elif each_taxid in changed_taxids and changed_taxids[each_taxid] in taxid_to_name_dict:
            all_taxons_print.append(":".join(
                [changed_taxids[each_taxid], taxid_to_name_dict[changed_taxids[each_taxid]].replace(":", "_").replace(",", "_").replace("\'", "")]))
        else:
            #print(f"Still cannot find {each_taxid}")
            taxid_to_name_dict[each_taxid] = "Unknown"
            all_missing_taxons_print.append(each_taxid)
    for i in range(0, len(all_taxons_print), 5000):
        print("\n".join(all_taxons_print[i: i + 5000]), file=open("%sTaxon_to_Species.txt" % (itol_output), "a+"))
    for i in range(0, len(all_missing_taxons_print), 5000):
        print("\n".join(all_missing_taxons_print[i: i + 5000]), file=open("%sMissing_Taxons.txt" % (itol_output), "a+"))

    #all_missing_taxons_print
    for i in range (0, len(itol_output_generic), 5000):
        print("\n".join(itol_output_generic[i:i+5000]), file=open(generic_filename, "a+"))
    for i in range (0, len(itol_output_readytogo), 5000):
        print("\n".join(itol_output_readytogo[i:i+5000]), file=open(specific_filename, "a+"))

    app.after(500, update_status, args=["Creating Histograms", "Yellow"])
    #time.sleep(30)
    try:
        app.after(1000, histo_maker, args=[output_dir, TM_thresolh, CxxH_threshold, histogram_data, TM_thresolh_upper,
                                           CxxH_threshold_upper, scatterplot_data])
    except Exception as e:
        app.after(1550, update_status, args=["ERROR creating Histograms", "Red"])


def last_bit():
    global Sequences, Haem_Weight_Dict, operon_type_dict, print_dict, count_dict, top_hits, TM_thresolh, CxxH_threshold,\
        output_file_name,TM_thresolh_upper, CxxH_threshold_upper, Haem_Weight_Ratio_To_Haem_Dict, scatterplot_data
    print(f"{datetime.now()}: Started making operon images.")
    histogram_data = list()
    mini_print = list()
    if Sequences:

        # triple_count_dict = dict()
        # operon_type_dict = dict()
        counter = 0
        #print(Haem_Weight_Dict)
        for key in Haem_Weight_Dict:

            if type(Haem_Weight_Dict[key]) == str:
                print(key, Haem_Weight_Dict[key], sep=" , ")
        listofTuples = sorted(Haem_Weight_Dict.items(), reverse=True, key=lambda x: x[1])

        # Iterate over the sorted sequence
        for elem in listofTuples:
            if float(upper_limit.value) >= float(elem[1]) >= float(lower_limit.value):
                #print(elem)
                key=elem[0]
                #print(elem[0] , " ::" , elem[1] )
                #while len(Haem_Weight_Dict) > 0:
                #highest_value = count_dict[max(Haem_Weight_Dict, key=Haem_Weight_Dict.get)]
                #key = max(Haem_Weight_Dict, key=Haem_Weight_Dict.get)
                if operon_type_dict[key] != "":
                    #print(operon_type_dict[key])
                    histogram_data.append("%s,%s" % (operon_type_dict[key], print_dict[key])) # eg "MtrAB,Helicobacter(1055530),NP_208019,NC_017376,WP_000959603(1254040_1254726)-0-0-25.37555000000002...."
                    try:
                        os.mkdir("%sTM%s_%s_Haem%s_%s\\Figures\\" % (output_dir, TM_Threshold.value, TM_Threshold_upper.value, Cxxh_Threshold.value, Cxxh_Threshold_upper.value))
                    except:
                        pass
                    image_dir2 = "%sTM%s_%s_Haem%s_%s\\Figures\\" % (output_dir, TM_Threshold.value, TM_Threshold_upper.value, Cxxh_Threshold.value, Cxxh_Threshold_upper.value)
                    for op_type in operon_type_dict[key].split(";"):
                        mini_print.append("%s,%s" % (op_type, print_dict[key]))
                        # print(op_type, print_dict[key], sep=",", file=open(str(output_dir) + str(Output_File.value), "a+"))
                        if counter < top_hits:
                            input_line = print_dict[key].split(",")
                            final_line = input_line[:3]
                            for each in range(3, len(input_line)):
                                final_line.append(input_line[each].split("(")[0] + input_line[each].split(")")[1])
                                final_line.append(input_line[each].split("(")[1].split(")")[0])
                            final_line = (",".join(final_line)).replace("(", ",").replace(")", "")
                            try:
                                Spoofed_genbank(final_line, op_type, image_dir2)
                            except Exception as e:
                                print("Error making image2")
                                print(e)
                        counter += 1
    else:
        highest_value = count_dict[max(count_dict, key=count_dict.get)]
        counter = 0
        #print(highest_value)
        for scores in range(highest_value + 1, -1, -1):
                ## nice and separated for easy reading
            for key in count_dict:
                if count_dict[key] == scores:  # check best operon first
                    #print(print_dict[key], file = output_Operons)  # print here, if you want ALL operons printed. or remove next 'if' line
                    if operon_type_dict[key] != "":
                        histogram_data.append("%s,%s" % (operon_type_dict[key], print_dict[key]))

                        try:
                            os.mkdir("%sTM%s_%s_Haem%s_%s\\Figures\\" % (output_dir, TM_Threshold.value,
                                                                        TM_Threshold_upper.value, Cxxh_Threshold.value,
                                                                        Cxxh_Threshold_upper.value))
                        except:
                            pass
                        image_dir2 = "%sTM%s_%s_Haem%s_%s\\Figures\\" % (output_dir, TM_Threshold.value,
                                                                        TM_Threshold_upper.value, Cxxh_Threshold.value,
                                                                        Cxxh_Threshold_upper.value)
                        for op_type in operon_type_dict[key].split(";"):
                            mini_print.append("%s,%s" % (op_type, print_dict[key]))
                            # print(op_type, print_dict[key], sep=",", file=open(str(output_dir) + str(Output_File.value), "a+"))
                            if counter < top_hits:
                                input_line = print_dict[key].split(",")
                                final_line = input_line[:3]
                                for each in range(3, len(input_line)):
                                    final_line.append(input_line[each].split("(")[0] + input_line[each].split(")")[1])
                                    final_line.append(input_line[each].split("(")[1].split(")")[0])
                                final_line = (",".join(final_line)).replace("(", ",").replace(")", "")
                                try:
                                    Spoofed_genbank(final_line, op_type, image_dir2)
                                except Exception as e:
                                    print("Error making image")
                                    print(e)
                            counter += 1
    for lines in range(0, len(mini_print), 5000):
        print("\n".join(mini_print[lines: lines + 5000]), sep="", file=open(output_file_name, "a+"))
    app.after(500, update_status, args=["Creating iTOL files...", "Yellow"])

    try:
        app.after(1000, make_itol_files, args=[output_dir, TM_thresolh, CxxH_threshold, histogram_data, TM_thresolh_upper,
                                           CxxH_threshold_upper, scatterplot_data, operon_type_dict])
    except Exception as e:
        app.after(1550, update_status, args=["ERROR creating iTOL", "Red"])
    del mini_print, operon_type_dict, print_dict, count_dict, Haem_Weight_Dict


def input_operons():
    app.after(100, update_status, args=["Looking for files...", "Yellow"])
    app.after(300, file_finder)


def output_image():
    app.after(100, update_status, args=["Looking for image output location...", "Yellow"])
    app.after(300, find_output_directory)


def find_output_directory():    # change where you output files to. Default is within this directory
    global output_directory_images
    global currdir
    tempdir = filedialog.askdirectory(initialdir=currdir, title='Please select a directory')
    # if len(tempdir) > 0:
    #    print("You chose " + str(tempdir).split("/")[-1])
    output_directory_images = tempdir + "\\"
    app.after(100, update_status, args=["Image output location updated.", "Green"])


def create_image():
    app.after(100, update_status, args=["Outputting image...", "Orange"])
    app.after(300, image_print, args=[operon_textified.value])

def clear_op():
    operon_textified.clear()


def image_print(operon_line):
    input_line = operon_line.replace(" ", "").replace("\n", "").replace("\t", ",").split(",")
    op_type = input_line[0]
    del input_line[0]
    final_line = input_line[:3]
    for each in range(3, len(input_line)):
        final_line.append(input_line[each].split("(")[0] + input_line[each].split(")")[1])
        final_line.append(input_line[each].split("(")[1].split(")")[0])
    final_line = (",".join(final_line)).replace("(", ",").replace(")", "")
    try:

        app.after(100, Spoofed_genbank, args=[final_line, op_type, output_directory_images])
        app.after(200, update_status, args=["Image output completed.", "Green"])
    except:
        app.after(200, update_status, args=["ERROR making image.", "Red"])


def file_finder():    #change where you get your sample files from. Default is within this directory
    global currdir
    global input_csv
    global cxxh_file
    global TM_File
    cxxh_file = ""
    TM_File = ""
    input_csv = ""

    filez = filedialog.askopenfilenames(initialdir=currdir, title='Select predicted operon file, Cxxh file and TM file')
    for i in range(0, len(filez)):
        # filenamesfortext.append(filez[i].split("/")[-1][:-4])  # JUST save file names (without extensions)
        # filenames.append(filez[i])  # save file names with extension and directory, for ease of searching later...redundant currenty
        BatchFile = filez[i]

        if "operons" in BatchFile.lower():
            input_csv = BatchFile
        elif "cx" in BatchFile.lower():
            cxxh_file = BatchFile
        elif "tm" in BatchFile.lower():
            TM_File = BatchFile
    app.after(100, update_status, args=["Files loaded", "Green"])


client_ref_text = Text(app, text="Output file name (will create CSV)", color="black", grid=[0, 0])  # [column, row]
Output_File = TextBox(app, text="Default", grid=[1, 0, 2, 1], width="fill")

#hair_bleached_text = Text(app, text="Proximate Proteins only", color="black", grid=[0, 1])  # [column, row]
#Proximate_Only = TextBox(app, text="No", grid = [1, 1,2,1], width="fill")
#previously_tested_text = Text(app, text="Dual Proteins only", color="black", grid=[0, 2])  # [column, row]
previously_tested_text = Text(app, text="Operon Type", color="black", grid=[0, 1])  # [column, row]
OperonType = ListBox(
    app, grid=[1, 1, 2, 1],
    width=160,
    height=70,
    items=["Any", "Cyc2", "MtrAB", "MtrCAB", "custom"],
    selected=["Any"], ##one to select by default
    scrollbar=True,
    multiselect=True)

#Dual_Proteins = TextBox(app, text="No", grid=[1, 2, 2, 1], width= "fill")
lab_ref_text = Text(app, text="TM Threshold", color="black", grid=[0, 2])  # [column, row]
TM_Threshold = TextBox(app, text="1", grid=[1, 2])
TM_Threshold_upper = TextBox(app, text="100", grid = [2, 2])
disclosed_medications_text = Text(app, text="Cxxh Threshold", color="black", grid=[0, 3])  # [column, row]
Cxxh_Threshold = TextBox(app, text="1", grid = [1,3])
Cxxh_Threshold_upper = TextBox(app, text="100", grid = [2,3])
top = Text(app, text="Top hits to print as image", color="black", grid=[0, 4])  # [column, row]
Top_Hits = TextBox(app, text="10", grid = [1,4,2,1], width="fill")
File_chooser = Text(app, text="Select file type (cvg/png/tig/jpg/gif/...)", color="black", grid=[0, 5])  # [column, row]
Filetype = TextBox(app, text="svg", grid = [1,5,2,1], width="fill")
Filesize_chooser = Text(app, text="Select file size", color="black", grid=[0, 6])  # [column, row]
Filesize = TextBox(app, text="25", grid = [1,6,2,1], width="fill")
lower_limit_text = Text(app, text="(kDa/Haem) ratio", color="black", grid=[0, 7])  # [column, row]
lower_limit = TextBox(app, text="0", grid=[1,7])
upper_limit = TextBox(app, text="10000", grid=[2,7])
itol_many_genome_files = ListBox(
    app, grid=[0, 9, 2, 1],
    width=160,
    height=35,
    items=["One genome file (iTOL)", "All genome files (iTOL)"],
    selected="All genome files (iTOL)", ##one to select by default
    scrollbar=True)
customs_value = TextBox(window, text="", grid=[0,0,3,1], width="fill")
button3 = PushButton(window, text="Add Operon", grid=[4, 0], command=add_customs)
upper_limit_text = Text(window, text="Operon of form \'H-P-F-O\'\nH=Haem\nP=Porin\nF=Fusion\nO=any", color="black", grid=[0, 1])  # [column, row]

hm_text = Text(app, text="HeatMap Colour Scale (below)", color="black", grid=[3, 0])  # [column, row]
operon_choice = ListBox(
    app, grid=[3, 1, 1, 9],
    width=160,
    height=250,
    items=['Blues','Greys', 'Purples', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
           'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
           'Pink', 'Spring', 'Summer', 'Autumn', 'Winter', 'cool', 'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'],
    selected="Blues", ##one to select by default
    scrollbar=True)
def fix_itol():
    global currdir
    itol_and_phy = []  # save directories for when you load them with other button
    app.after(200, update_status, args=["Loading files...", "White"])
    filez = filedialog.askopenfilenames(initialdir=currdir, title='Choose sample files')
    for i in range(0, len(filez)):
        # filenamesfortext.append(filez[i].split("/")[-1][:-4])  # JUST save file names (without extensions)
        itol_and_phy.append(
            filez[i])  # save file names with extension and directory, for ease of searching later...redundant currenty
        # pdf_file = filez[i]
    all_taxons = []
    for filename in filez:
        seq_file=open(filename, "r")
        for line in seq_file:
            all_taxons.append(line.split(",")[0].split("(")[-1].strip(")"))
    all_taxons_set =list(set(all_taxons))
    all_taxons_print = []
    all_taxids_dict = {}
    for line in open(currdir+"taxdmp\\names.dmp", "r"):
        data = line.split("|")
        all_taxids_dict[data[0].strip()] = data[1]
    for each_taxid in all_taxons_set:
        if each_taxid in all_taxids_dict:
            all_taxons_print.append(":".join(each_taxid, all_taxids_dict[each_taxid]))
        else:
            print("MISSING TAXID?: " + each_taxid)
    for i in range(0, len(all_taxons_print), 5000):
        print("\n".join(all_taxons_print[i: i+5000]), file=open("Taxid_to_Species.txt", "a+"))
    # info("Data loaded", "All data has loaded successfully.")
    app.after(500, update_status, args=["iTOL files fixed.", "Green"])


operon_textified = TextBox(window2, text="", grid=[0,0,3,1], width="fill")
run_operon = PushButton(window2, text="Create image", grid=[4, 0], command=create_image)
clear_operon = PushButton(window2, text="Clear", grid=[4, 1], command=clear_op)
operon_words = Text(window2, text="Tabs = commas, spaces ignored.\nExample format\n(replace newline with commas):\nMtrAB\nGeobacter(1340425)\nNP_953690\nNZ_CP014963\nProt1*(4_100)-12-0-33.31\nProt2(150_210)-0-20-47.07\nProt3(220_400)-5-0-18.03\nProt4(420_469)-6-0-24.72", color="black", grid=[0, 1])  # [column, row]


hide_query_duplicate = ListBox(
    app, grid=[1, 9, 2, 1],
    width=160,
    height=35,
    items=["Hide duplicate query ID", "Do not hide"],
    selected="Do not hide", ##one to select by default
    scrollbar=True)

menubar = MenuBar(app,
                  toplevel=["File"],
                  options=[
                      [["Input Operon, TM and Cxxh files", input_operons]]])
#menubar = MenuBar(app, toplevel=["File", "iTOL"], options=[ [["Input Operon, TM and Cxxh files", input_operons]], [["Input Phy and Taxon files", fix_itol]]])

menubar_w2 = MenuBar(window2,
                  toplevel=["File"],
                  options=[
                      [["Output location", output_image]]])

button1 = PushButton(app, text="Custom operon search", grid=[0, 10], command=windowed)
button_image = PushButton(app, text="Create image", grid=[1, 10], command=windowed2)
button2 = PushButton(app, text="Run app", grid=[2, 10], command=Classifyer)
status = Text(app, text="Status: Ready", grid=[3, 10])

app.display()
