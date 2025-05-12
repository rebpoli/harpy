#!/usr/bin/env -S python

import pandas as pd
import matplotlib.pyplot as plt
import os, time
from matplotlib.animation import FuncAnimation

fn = "py/paper.mplstyle"
if os.path.isfile(fn) :
    plt.style.use('default')   ## reset!
    plt.style.use(fn)

def plot_csv( fn ) :

    plt.show()

class CSVPlotter:
    def __init__(self, csv_file, refresh_interval=1000):
        self.csv_file = csv_file
        self.refresh_interval = refresh_interval
        self.last_modified = 0

        # Set up the plot
        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        plt.title('Live CSV Data Visualization')

        # Initial plot
        self.update(0)

        # Set up animation
        self.ani = FuncAnimation(
            self.fig, self.update,
            interval=self.refresh_interval,
            save_count=10, blit=False
        )

    #
    #
    def do_plot( self ) :
        df = pd.read_csv(self.csv_file, sep="\t")
        df = df[df["Time(s)"]>0]

        # Index(['Timestep', 'Time(s)', 'Time(day)', 'Var', 'X', 'Y', 'Z', 'Value'], dtype='object')
        ax = self.ax

        sxxdf = df[ df.Var == "STOTXX" ]
        sxxdf = sxxdf[ sxxdf.X > -0.02 ]
        sxxdf = sxxdf[ sxxdf.X <  0 ]

        ##
        ax.plot( sxxdf["Time(day)"], sxxdf["Value"], marker='x', markersize=2)
        ##

        ax.set_xscale('log')
        ax.set_xlabel("Time (days)")
        ax.set_ylabel(r"Stress at the sphere center ($\sigma_{xx}$) - MPa")

    #
    #
    def update(self, frame):
        if not os.path.exists( self.csv_file ) : return;

        current_mtime = os.path.getmtime(self.csv_file)

        if current_mtime > self.last_modified:
            self.last_modified = current_mtime
            print(f"CSV updated at {time.strftime('%H:%M:%S')} - refreshing plot")

            # Clear the current plot
            self.ax.clear()

            # Read and plot the data
            try:
                self.do_plot()
            except Exception as e:
                print(f"Error updating plot: {e}")

        return self.ax,

    #
    #
    def show(self):
        """Display the plot window"""
        plt.tight_layout()
        plt.show()

#
#
#
#
csv_file = "run/csv/lin_0.csv"
plotter = CSVPlotter(csv_file, refresh_interval=1000)  # Check every 1000ms (1 second)
print(f"Monitoring {csv_file}. Close the plot window to stop.")
plotter.show()
