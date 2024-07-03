from tkinter import *
from tkinter import scrolledtext

def show_gui(protein_seq, zf_list, predictions, threshold=0.5):
    # Function to create the Tkinter GUI for displaying protein sequence and highlighting domains
    root = Tk()
    root.geometry("900x600")
    root.title("Protein Sequence Viewer")
    root.resizable(0, 0)

    label = Label(root, text="Protein Sequence Viewer", font=('Times New Roman', 16))
    label.pack(pady=10)

    text_area = scrolledtext.ScrolledText(root, width=100, height=30, wrap=WORD, font=("Courier", 12))
    text_area.pack(expand=True, fill="both", padx=20, pady=10)

    # Split the protein sequence into lines and add a line space above each line
    lines = protein_seq.split('\n')
    spaced_protein_seq = '\n\n'.join(lines)

    text_area.insert(END, spaced_protein_seq)

    # Highlight ZF domains in different colors and display predictions
    for zf_seq, prediction in zip(zf_list, predictions):
        start = spaced_protein_seq.find(zf_seq)
        end = start + len(zf_seq)

        while start != -1:
            text_area.tag_add("zf", f"1.{start}", f"1.{end}")
            text_area.tag_config("zf", background="yellow", foreground="black")

            # Display prediction next to the zinc finger domain
            text_area.insert(f"1.{end}", f"({prediction}) ", "prediction")
            text_area.tag_config("prediction", foreground="orange")

            start = spaced_protein_seq.find(zf_seq, start + 1)
            end = start + len(zf_seq)

    root.mainloop()

def get_input_and_display():
    # Function to get input from user and display the GUI
    protein_seq = input_seq.get("1.0", "end-1c")
    threshold = float(input_threshold.get()) if input_threshold.get().strip() else 0.5  # Default threshold is 0.5 if not provided

    # Example list of extracted zinc finger domains and predictions
    extracted_zf_list = [
        'GKFFSWRSNLTR', 'GKSFSRSSHLIG', 'GKSFSWFSHLVT', 'GKSFVHSSRLIR',
        'GKSFRQSTHLIL', 'GKSYSQRSHLVV', 'GKCFSRSSHLYS', 'GKSFSQSSALIV',
        'GKAFIRKNDLIK', 'GIIFSQNSPFIV', 'GTALVNTSNLIG'
    ]

    predictions = [
        '0.8', '0.7', '0.6', '0.61',
        'Low', 'High', 'Low', 'Medium',
        'Low', 'High', 'Medium'
    ]

    # Clear the input fields
    input_seq.delete("1.0", "end")
    input_threshold.delete(0, "end")

    # Call function to create GUI with protein sequence, highlighted domains, and predictions
    show_gui(protein_seq, extracted_zf_list, predictions, threshold)

# Main Tkinter window setup
win = Tk()
win.geometry("700x400")
win.title("Protein Sequence Viewer Setup")
win.resizable(0, 0)

# Label and entry for protein sequence input
label_seq = Label(win, text="Enter protein sequence:")
label_seq.pack(pady=10)
input_seq = Text(win, height=5, width=60, wrap=WORD)
input_seq.pack(pady=10)

# Label for threshold input (optional)
label_threshold = Label(win, text="Enter threshold (optional):")
label_threshold.pack()

# Entry for threshold input
input_threshold = Entry(win, width=30)
input_threshold.pack(pady=10)

# Button to submit inputs and display GUI
submit_button = Button(win, text="Submit", command=get_input_and_display)
submit_button.pack()

win.mainloop()


# protein_seq = "MDAKSLTAWSRTLVTFKDVFVDFTREEWKLLDTAQQIVYRNVMLENYKNLVSLGYQLTKPDVILRLEKGEEPWLVEREIHQETHPDSETAFEIKSSVSSRSIFKDKQSCDIKMEGMARNDLWYLSLEEVWKCRDQLDKYQENPERHLRQVAFTQKKVLTQERVSESGKYGGNCLLPAQLVLREYFHKRDSHTKSLKHDLVLNGHQDSCASNSNECGQTFCQNIHLIQFARTHTGDKSYKCPDNDNSLTHGSSLGISKGIHREKPYECKECGKFFSWRSNLTRHQLIHTGEKPYECKECGKSFSRSSHLIGHQKTHTGEEPYECKECGKSFSWFSHLVTHQRTHTGDKLYTCNQCGKSFVHSSRLIRHQRTHTGEKPYECPECGKSFRQSTHLILHQRTHVRVRPYECNECGKSYSQRSHLVVHHRIHTGLKPFECKDCGKCFSRSSHLYSHQRTHTGEKPYECHDCGKSFSQSSALIVHQRIHTGEKPYECCQCGKAFIRKNDLIKHQRIHVGEETYKCNQCGIIFSQNSPFIVHQIAHTGEQFLTCNQCGTALVNTSNLIGYQTNHIRENAY"
