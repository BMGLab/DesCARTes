import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import io
import os

# The data string from the table
csv_data = """Sayı,Hasta,Tümör Bloğu MUC1,Tümör Bloğu CLAUDIN-4,Sağlıklı eşlenik doku bloğu MUC1,Sağlıklı eşlenik doku bloğu CLAUDIN-4
1,2043-24,DDPOL,2+,Negatif,Negatif
2,5337-24,DEPOL,3+,Negatif,Negatif
3,6372-24,YDPOL,1+,Negatif,Negatif
4,10274-24,YDPOL,3+,Negatif,Negatif
5,11172-24,DDPOL,2+,Negatif,Negatif
6,11173-24,DEPOL,1+,Negatif,Negatif
7,11700-24,YDPOL,1+,Negatif,Negatif
8,11702-24,DEPOL,2+,Negatif,Negatif
9,14553-24,DDPOL,2+,Negatif,Negatif
10,16102-24,1-,1+,Negatif,Negatif
11,19405-24,DEPOL,2+,Negatif,Negatif
12,19811-24,DEPOL,1+,Negatif,Negatif
13,20680-24,YDPOL,3+,Negatif,Negatif
14,20887-24,DEPOL,2+,Negatif,Negatif
15,20989-24,YDPOL,3+,Negatif,Negatif
16,24184-24,DEPOL,3+,Negatif,Negatif
17,24834-24,1-,1+,Negatif,Negatif
18,24999-24,DEPOL,3+,Negatif,Negatif
19,25241-24,DEPOL,2+,Negatif,Negatif
20,27685-24,1-,1+,Negatif,Negatif
21,28255-24,YDPOL,2+,Negatif,Negatif
22,28593-24,DDPOL,1+,Negatif,Negatif
23,28784-24,YDPOL,3+,Negatif,Negatif
24,29369-24,DEPOL,2+,Negatif,Negatif
25,29454-24,DEPOL,3+,Negatif,Negatif
26,29457-24,YDPOL,1+,Negatif,Negatif
27,32481-24,DEPOL,2+,Negatif,Negatif
28,32648-24,DDPOL,2+,Negatif,Negatif
29,32815-24,DDPOL,2+,Negatif,Negatif
30,33101-24,1-,1+,Negatif,Negatif
31,33472-24,YDPOL,2+,Negatif,Negatif
32,33687-24,DDPOL,3+,Negatif,Negatif
33,33915-24,YDPOL,1+,Negatif,Negatif
34,37554-24,YDPOL,3+,Negatif,Negatif
35,37876-24,YDPOL,1+,Negatif,Negatif
36,37993-24,DEPOL,2+,Negatif,Negatif
37,38176-24,DDPOL,2+,Negatif,Negatif
38,46048-24,1-,1+,Negatif,Negatif
39,46290-24,DEPOL,1+,Negatif,Negatif
40,46719-24,1-,1+,Negatif,Negatif
41,47077-24,YDPOL,3+,Negatif,Negatif
42,47398-24,DEPOL,3+,Negatif,Negatif
43,47542-24,DEPOL,1+,Negatif,Negatif
44,47705-24,YDPOL,3+,Negatif,Negatif
45,48287-24,YDPOL,3+,Negatif,Negatif"""

# Load the data into a DataFrame
df = pd.read_csv(io.StringIO(csv_data))

# Map the raw expression scores to match the specific color scheme keys
mapping = {
    '1+': 'Low', 
    '2+': 'Medium', 
    '3+': 'High',
    '1-': 'Not detected',
    'Negatif': 'Not detected'
}
df['CLAUDIN-4 Category'] = df['Tümör Bloğu CLAUDIN-4'].map(mapping)

# Define your requested color scheme
custom_colors = {
    "High": "#E69F00",        # orange
    "Medium": "#56B4E9",      # sky blue
    "Low": "#009E73",         # green
    "Not detected": "#ECECEC" # light gray
}

# Ask the user to define an output path
output_dir = input("Please enter the directory path where you want to save the plots (e.g., C:/Users/Name/Documents or ./outputs): ")

# Ensure the directory exists, create it if it doesn't
os.makedirs(output_dir, exist_ok=True)

# Set a standard sans-serif font for matplotlib to closely match the image
plt.rcParams['font.family'] = 'sans-serif'

# Create the plot
plt.figure(figsize=(8, 6))

# Plotting a countplot and assigning it to 'ax'
ax = sns.countplot(
    data=df, 
    x='CLAUDIN-4 Category', 
    order=['Not detected', 'Low', 'Medium', 'High'], # Keeps the x-axis logical
    palette=custom_colors
)

# Add 'n=' labels to the top of each bar
for container in ax.containers:
    ax.bar_label(
        container, 
        fmt='n=%d',             # Prepend n= before the integer count
        padding=4, 
        fontsize=12,            # Adjusted size to match the visual proportion
        color='#595959',        # Dark gray color similar to the screenshot
        fontweight='normal'     # Regular weight
    )

# Format the titles and labels
plt.title('CLAUDIN-4 Expression Level Distribution in Tumor Samples', fontsize=14, pad=15)
plt.xlabel('Expression Level', fontsize=12)
plt.ylabel('Number of Patients', fontsize=12)

# Adjust y-axis limit slightly so the new labels don't get cut off at the top
max_count = df['CLAUDIN-4 Category'].value_counts().max()
plt.ylim(0, max_count * 1.15) 

# Adjust layout to make sure labels fit
plt.tight_layout()

# Construct file paths for both PNG and PDF
png_path = os.path.join(output_dir, 'claudin4_IHC_bar_chart.png')
pdf_path = os.path.join(output_dir, 'claudin4_IHC_bar_chart.pdf')

# Save the plots
plt.savefig(png_path, format='png', dpi=300) # Added high DPI for better quality
plt.savefig(pdf_path, format='pdf')

print(f"\nPlots successfully saved to:\n- {png_path}\n- {pdf_path}")

# Optionally show the plot on screen
plt.show()
