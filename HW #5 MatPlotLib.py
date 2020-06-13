#!/usr/bin/env python
# coding: utf-8

# In[89]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np

# Study data files
mouse_metadata_path = r"C:\Users\jimkn\Downloads\Homework_05-Matplotlib_Instructions_Pymaceuticals_data_Mouse_metadata.csv"
study_results_path = r"C:\Users\jimkn\Downloads\Homework_05-Matplotlib_Instructions_Pymaceuticals_data_Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Merge data sets
merged_df = pd.merge(mouse_metadata, study_results, on="Mouse ID",how="left")
merged_df.head()


# In[ ]:





# In[33]:


# Generate a summary statistics table of mean, median, variance,
#standard deviation, and SEM of the tumor volume for each regimen

means = merged_df.groupby('Drug Regimen').mean()['Tumor Volume (mm3)']
medians = merged_df.groupby('Drug Regimen').median()['Tumor Volume (mm3)']
variances = merged_df.groupby('Drug Regimen').var()['Tumor Volume (mm3)']
standards = merged_df.groupby('Drug Regimen').std()['Tumor Volume (mm3)']
sems = merged_df.groupby('Drug Regimen').sem()['Tumor Volume (mm3)']

newtable = pd.DataFrame(means)
newtable2 = newtable.rename(columns={"Tumor Volume (mm3)": "Mean"})

newtable2["Median"] = medians
newtable2["Variance"] = variances
newtable2["std"] = standards
newtable2["sem"] = sems

newtable2


# In[10]:


# Bar Plots


# In[34]:


# Generate a bar plot showing number of data points for each treatment regimen using pandas

data = merged_df.groupby('Drug Regimen').count()['Tumor Volume (mm3)']
treatment_bar = pd.DataFrame(data)

bar_graph = treatment_bar.plot.bar(legend=True)
bar_graph
plt.ylabel("Number of Data Points")
plt.title("Data Points Per Drug Treatment Regimen")
plt.savefig('barplot1')


# In[35]:


# Generate a bar plot showing number of data points for each treatment regimen using pyplot

x_axis = np.arange(len(data))

tick_locations = [x for x in x_axis]

plt.bar(x_axis, treatment_bar['Tumor Volume (mm3)'], alpha=0.75, align="center")
plt.xticks(tick_locations, newtry['Drug Regimen'],rotation="vertical")

plt.xlim(-0.75, len(data)-.25)
plt.ylim(0, 250)

plt.title("Data Points Per Drug Treatment Regimen")
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Data Points")

plt.savefig('barplot2')
plt.show()


# In[26]:


# Pie Plots


# In[53]:


# Create Variables for Pie Plots
bygender = mouse_metadata.groupby("Sex").count()

labels = [mouse_metadata['Sex'].unique()]
newlist = list(bygender.index)

sizes = [bygender["Mouse ID"]]


# In[54]:


# Generate a pie plot showing the distribution of female versus male mice using pandas

colors = ["orange", "blue"]

plt.pie(sizes, labels=newlist, colors=colors, autopct="%1.1f%%")
plt.title('Male vs Female Mouse Population')
plt.ylabel('Sex Percentage')

plt.show()


# In[60]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot

secondpie = bygender["Mouse ID"].plot.pie(y='sizes', autopct="%1.1f%%")
plt.title('Male vs Female Mouse Population')
plt.ylabel('Sex Percentage')

plt.show()


# In[61]:


# Quartiles, Outliers and Boxplots


# In[64]:


# Calculate the final tumor volume of each mouse across four of the most promising treatment regimens. Calculate the IQR and quantitatively determine if there are any potential outliers
#Extract the top 4 regimenes from the data frame in order to perform IQR test on each
best_regimes = merged_df[merged_df["Drug Regimen"].isin(["Capomulin", "Ramicane", "Infubinol", "Ceftamin"])]
best_regimes = best_regimes.sort_values(["Timepoint"], ascending=True)
best_regimes

best_regimes_data = best_regimes[["Drug Regimen", "Mouse ID", "Timepoint", "Tumor Volume (mm3)"]]

best_regimes_data


# In[67]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest

# Find last tumor treatment
last_treatment = best_regimes_data.groupby(['Drug Regimen', 'Mouse ID']).last()['Tumor Volume (mm3)']
last_treatment


# In[69]:


last_treatment_df = last_treatment.to_frame()
last_treatment_df


# In[76]:


# Create label for boxplot
top_4 = ['Capomulin', 'Ramicane', 'Infubinol','Ceftamin']

# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
final_df = last_treatment_df.reset_index()
tumor_lists = final_df.groupby('Drug Regimen')['Tumor Volume (mm3)'].apply(list)
tumor_list_df = pd.DataFrame(tumor_lists)
tumor_list_df = tumor_list_df.reindex(top_4)
tumor_volumes = [volumes for volumes in tumor_list_df['Tumor Volume (mm3)']]
plt.boxplot(tumor_volumes, labels=top_4)
plt.ylim(10, 80)
plt.show()


# In[77]:


# Line and Scatter Plots


# In[88]:


# Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin
# Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin
time_vs_tumor = merged_df[merged_df["Mouse ID"].isin(["j119"])]
time_vs_tumor

time_vs_tumor_data = time_vs_tumor[["Mouse ID", "Timepoint", "Tumor Volume (mm3)"]]
time_vs_tumor_data

line_plot_df = time_vs_tumor_data.reset_index()
line_plot_df

line_plot_final = line_plot_df[["Mouse ID", "Timepoint", "Tumor Volume (mm3)"]]
line_plot_final

lines = line_plot_final.plot.line()


# In[82]:


# Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin regimen
capomulin_scatter = merged_df[merged_df["Drug Regimen"].isin(["Capomulin"])]

capomulin_scatter_df = best_regimes[["Mouse ID","Weight (g)", "Tumor Volume (mm3)"]]

capomulin_sorted = capomulin_scatter_df.sort_values(["Weight (g)"], ascending=True)

capomulin_scatter_plot = capomulin_sorted.reset_index()

capomulin_grouped_weight = capomulin_scatter_plot.groupby("Weight (g)")["Tumor Volume (mm3)"].mean()

capo_grouped_plot = pd.DataFrame(capomulin_grouped_weight).reset_index()

capomulin_scatter = capo_grouped_plot.plot(kind='scatter', x='Weight (g)', y='Tumor Volume (mm3)', grid = True, figsize= (8,8))


# In[83]:


# Calculate the correlation coefficient and linear regression model for mouse weight and average tumor volume for the Capomulin regimen


# In[86]:


x_values = capo_grouped_plot["Weight (g)"]
y_values = capo_grouped_plot["Tumor Volume (mm3)"]
linregress = (x_values, y_values)
(slope, intercept, rvalue, pvalue, stderr) = linregress(x_values, y_values)
regress_values = x_values * slope + intercept
line_eq = "y =" + str(round(slope,2)) + "x + " + str(round(intercept,2))
plt.scatter(x_values, y_values)
plt.plot(x_values,regress_values,"r-")
plt.annotate(line_eq,(6,10),fontsize=10,color="red")
plt.xlabel("Weight")
plt.ylabel("Tumor Volume")
plt.title("Weight Vs. Avg Tumor Vol")
plt.show()


# In[ ]:




