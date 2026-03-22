import os
import re
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Import BioPython modules for fetching real gene data.
from Bio import Entrez, SeqIO
Entrez.email = "omar.m.abdalaal@gmail.com"  # Replace with your email

########################################
# Fetch Real Gene Sequences from NCBI  #
########################################

# Define a dictionary of gene names along with their corresponding accession numbers.
# (These accession numbers are provided as examples. You may wish to verify or update them.)
gene_accessions = {
    "BDNF": "NM_170731.3",   # Brain-Derived Neurotrophic Factor
    "NRG1": "NM_013964.4",   # Neuregulin 1
    "DRD2": "NM_000795.4",   # Dopamine Receptor D2
    "GRIA1": "NM_000831.3"   # Glutamate Ionotropic Receptor AMPA type subunit 1
}

def fetch_gene_sequences():
    """
    For each gene in the gene_accessions dictionary, fetch its real sequence data
    from NCBI (using the Entrez API). Returns a dictionary compatible with our GUI.
    """
    genes = {}
    for gene, acc in gene_accessions.items():
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            genes[gene] = {
                "description": record.description,
                "sequence": str(record.seq.upper())
            }
            handle.close()
            print(f"Fetched data for gene: {gene}")
        except Exception as e:
            print(f"Error fetching gene {gene} (Accession: {acc}): {e}")
    return genes

#########################################
# Helper Functions: Motif Analysis, EEG #
#########################################

def search_motif(sequence, motif):
    """
    Searches for all occurrences of a motif in a given sequence.
    Returns:
        count: number of matches
        positions: list of starting indices (0-indexed) where the motif is found.
    The search is case-insensitive.
    """
    motif = motif.upper()
    sequence = sequence.upper()
    positions = [m.start() for m in re.finditer('(?={})'.format(re.escape(motif)), sequence)]
    count = len(positions)
    return count, positions

def simulate_eeg(motif_count, duration=5, fs=256):
    """
    Simulates an EEG signal.
    
    The EEG signal is modeled as:
      - A baseline alpha wave (10 Hz sine wave)
      - A gamma wave (40 Hz sine wave) whose amplitude is modulated by motif_count.
      - Added random noise.
    
    Parameters:
      motif_count: how many times the motif was found (affects gamma amplitude).
      duration: Duration of the signal in seconds.
      fs: Sampling frequency in Hz.
      
    Returns:
      t: time vector.
      eeg_signal: simulated EEG signal as a numpy array.
    """
    t = np.linspace(0, duration, int(fs * duration), endpoint=False)
    
    # Baseline alpha component: 10 Hz sine wave.
    alpha_wave = 50 * np.sin(2 * np.pi * 10 * t)
    
    # Gamma component: 40 Hz sine wave. Its amplitude increases with the motif count.
    gamma_amplitude = 1 + 0.5 * motif_count
    gamma_wave = gamma_amplitude * np.sin(2 * np.pi * 40 * t)
    
    # Add some Gaussian noise.
    noise = np.random.normal(0, 5, len(t))
    
    eeg_signal = alpha_wave + gamma_wave + noise
    return t, eeg_signal

#########################################
# Graphical User Interface (GUI) Design #
#########################################

class IntegratedEEGApp:
    def __init__(self, master):
        self.master = master
        master.title("Neural Gene Motif Finder & EEG Simulator")
        
        # Fetch gene data from NCBI.
        self.genes = fetch_gene_sequences()
        if not self.genes:
            messagebox.showerror("Data Error", "Unable to fetch gene sequences. Please check your internet connection or your API settings.")
            master.quit()
            return
        
        self.gene_names = list(self.genes.keys())
        
        # Create Frames
        self.top_frame = ttk.Frame(master)
        self.top_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)
        
        self.bottom_frame = ttk.Frame(master)
        self.bottom_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Gene Selection
        ttk.Label(self.top_frame, text="Select Gene: ").grid(row=0, column=0, sticky=tk.W)
        self.gene_var = tk.StringVar(value=self.gene_names[0])
        self.gene_dropdown = ttk.OptionMenu(self.top_frame, self.gene_var, self.gene_names[0], *self.gene_names, command=self.update_gene_info)
        self.gene_dropdown.grid(row=0, column=1, sticky=tk.W)
        
        # Gene Description Display
        self.gene_desc_label = ttk.Label(self.top_frame, text="Gene Description: ", wraplength=400)
        self.gene_desc_label.grid(row=1, column=0, columnspan=3, sticky=tk.W, pady=(5, 0))
        self.update_gene_info(self.gene_names[0])
        
        # Motif Entry
        ttk.Label(self.top_frame, text="Enter DNA Motif: ").grid(row=2, column=0, sticky=tk.W, pady=(5, 0))
        self.motif_entry = ttk.Entry(self.top_frame, width=20)
        self.motif_entry.grid(row=2, column=1, sticky=tk.W, pady=(5, 0))
        
        # Run Analysis Button
        self.run_button = ttk.Button(self.top_frame, text="Run Analysis", command=self.run_analysis)
        self.run_button.grid(row=2, column=2, padx=10, pady=(5, 0))
        
        # Results Text Area
        self.result_text = tk.Text(self.top_frame, height=4, width=60)
        self.result_text.grid(row=3, column=0, columnspan=3, pady=(10, 0))
        
        # Set up the Matplotlib figure for EEG simulation.
        self.fig, self.ax = plt.subplots(figsize=(6, 3))
        self.ax.set_title("Simulated EEG Signal")
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Amplitude")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.bottom_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def update_gene_info(self, selected_gene):
        """Update the description label based on the selected gene."""
        gene_info = self.genes.get(selected_gene, {})
        desc = gene_info.get("description", "No description available")
        self.gene_desc_label.config(text=f"Gene Description: {desc}")
    
    def run_analysis(self):
        """Perform motif search and simulate EEG signal; update GUI results."""
        selected_gene = self.gene_var.get()
        motif = self.motif_entry.get().strip()
        
        if not motif:
            messagebox.showerror("Input Error", "Please enter a DNA motif.")
            return
        
        gene_data = self.genes.get(selected_gene)
        if not gene_data:
            messagebox.showerror("Gene Error", "Selected gene not found.")
            return
        
        sequence = gene_data["sequence"]
        # Perform motif search.
        count, positions = search_motif(sequence, motif)
        
        # Update the result text area.
        result_message = (
            f"Selected Gene: {selected_gene}\n"
            f"Motif: {motif}\n"
            f"Motif Count: {count}\n"
            f"Positions: {positions if positions else 'None found'}\n"
            f"Sequence Length: {len(sequence)} nucleotides"
        )
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, result_message)
        
        # Simulate the EEG signal using the motif count.
        t, eeg_signal = simulate_eeg(count)
        
        # Update the EEG plot.
        self.ax.clear()
        self.ax.plot(t, eeg_signal, label="EEG Signal", color="blue")
        self.ax.set_title("Simulated EEG Signal")
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Amplitude")
        self.ax.legend()
        self.canvas.draw()

##########################
# Start the Tkinter App  #
##########################

if __name__ == "__main__":
    root = tk.Tk()
    app = IntegratedEEGApp(root)
    root.mainloop()