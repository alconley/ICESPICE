import polars as pl
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

tex_fonts = {
    # Use LaTeX to write all text
    # "text.usetex": True,
    "font.family": "serif",
    "font.serif" : ["CMR10"],
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 40,
    "font.size": 40,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 35,
    "xtick.labelsize": 35,
    "ytick.labelsize": 35
}


plt.rcParams.update(tex_fonts)
matplotlib.rcParams['axes.unicode_minus'] = False


df_withICESPICE = pl.read_parquet("./207Bi/207Bi_ICESPICE_f70mm_g30mm_run_*.parquet")
df_withoutICESPICE = pl.read_parquet("./207Bi/207Bi_noICESPICE_f9mm_g0mm_run_13.parquet")
# create a energy calibrated column of PIPS1000Energy with m=0.5395 and b=2.5229

df_withICESPICE = df_withICESPICE.with_columns([
    (pl.col("PIPS1000Energy") * 0.5395 + 2.5229).alias("PIPS1000EnergyCalibrated")
])

df_withoutICESPICE = df_withoutICESPICE.with_columns([
    (pl.col("PIPS1000Energy") * 0.5395 + 2.5229).alias("PIPS1000EnergyCalibrated")
])

fig, axs = plt.subplots(1, 1, figsize=(22, 6))
fig.subplots_adjust(left=0.075, right=0.99, top=0.99, bottom=0.2)

# scale of the hist of the data with out ICESPICE so the counts match the 975 peak
scale = 51700/65625

axs.hist(df_withoutICESPICE["PIPS1000EnergyCalibrated"], bins=1000, range=[200, 1200], histtype="step", color='#5CB8B2', weights=[scale]*len(df_withoutICESPICE), label="without ICESPICE", linewidth=2)

# plot the histogram of the calibrated energy column
axs.hist(df_withICESPICE["PIPS1000EnergyCalibrated"], bins=1000, range=[200, 1200], histtype="step", color="#A6192E", label="with ICESPICE", linewidth=2)


# label the lines
fs = 25
axs.text(481, 3300, r"570K", horizontalalignment='center', verticalalignment='bottom', fontsize=fs, rotation=90)
axs.text(553, 1500, r"570L", horizontalalignment='center', verticalalignment='bottom', fontsize=fs, rotation=90)
axs.text(565, 1300, r"570M", horizontalalignment='left', verticalalignment='bottom', fontsize=fs, rotation=90)
axs.text(975, 5150, r"1064K", horizontalalignment='center', verticalalignment='bottom', fontsize=fs, rotation=90)
axs.text(1047, 1400, r"1064L", horizontalalignment='center', verticalalignment='bottom', fontsize=fs, rotation=90)
axs.text(1059, 700, r"1064M", horizontalalignment='left', verticalalignment='bottom', fontsize=fs, rotation=90)

axs.set_ylim(0, 6500)

# set the labels
axs.set_xlabel(r"Energy [keV]")
axs.set_ylabel(r"Counts/keV")

# set the legend
axs.legend(loc='upper left', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')

# nice ticks
axs.minorticks_on()
axs.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=10)
axs.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=20)

# save the plot
# plt.savefig("./207Bi/207Bi_ICESPICE_spectrum.png", dpi=300)

# show the plot
plt.show()


# # Create the histograms using numpy.histogram with 1200 bins and range 0-1200
# hist_withICESPICE, bin_edges_with = np.histogram(
#     df_withICESPICE["PIPS1000EnergyCalibrated"], 
#     bins=1200, 
#     range=[0, 1200]
# )

# hist_withoutICESPICE, bin_edges_without = np.histogram(
#     df_withoutICESPICE["PIPS1000EnergyCalibrated"], 
#     bins=1200, 
#     range=[0, 1200], 
#     weights=[scale] * len(df_withoutICESPICE)
# )

# # Export the histogram data for ROOT macro or any other format
# # Save histogram bins and counts for with ICESPICE
# with open("207Bi_withICESPICE.txt", "w") as f:
#     for i in range(len(hist_withICESPICE)):
#         f.write(f"{bin_edges_with[i]:.2f} {hist_withICESPICE[i]:.2f}\n")

# # Save histogram bins and counts for without ICESPICE
# with open("207Bi_withoutICESPICE.txt", "w") as f:
#     for i in range(len(hist_withoutICESPICE)):
#         f.write(f"{bin_edges_without[i]:.2f} {hist_withoutICESPICE[i]:.2f}\n")
        