print "Initializing..."
import xlrd
import scipy as sp 
#import math
import numpy as np
import matplotlib.pyplot as plt
print "Initialization complete."

print "Begin calculation output: "
  
# Define the path to the .xls file(s) containing the conversion data.
# import the worksheet as a sheet object

data_worksheet = xlrd.open_workbook(filename=data_path[i]).sheet_by_index(0)
    
# Import conversion data from worksheet and store as scipy arrays
flow = sp.array(data_worksheet.col_values(0, start_rowx=4, end_rowx=None))

# Generate and annotate plots
plt.figure(1)
# Change labels, linestyle, and marker style here
plt.plot(flow0,HCeff0,'k*',label='NW with 20 nm Pt')
plt.xlabel("Flow Rate (sccm)") 
plt.ylabel("HC Conversion Efficiency (%)")
plt.title("HC Conversion Efficiency v. Flow Rate")
plt.ylim(0,100)
plt.legend()
plt.grid(True)
    
plt.savefig('Plots/deltaHC.pdf')
plt.savefig('Plots/deltaHC.png')
plt.show()

print "\nCalculation output completed."