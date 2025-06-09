import pandas as pd

#HTML template for the report
HTML_TEMPLATE = """
<html>
<head>
    <title>VCF ANALYSIS REPORT</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 20px;
            background-color: #f7f9fc;
            color: #333;
        }}
        h2 {{
            text-align: center;
            text-transform: uppercase;
            font-size: 2.5em;
            margin-bottom: 30px;
            color: #2c3e50;
            letter-spacing: 3px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            max-width: 100%;
            background-color: #ffffff;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px 15px;
            text-align: left;
        }}
        th {{
            background-color: #4a90e2;
            color: white;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        tr:nth-child(even) {{
            background-color: #f2f6fc;
        }}
        tr:hover {{
            background-color: #dbe9ff;
        }}
    </style>
</head>
<body>
    <h2>VCF ANALYSIS REPORT</h2>
    {table}
</body>
</html>
"""


#Function to generate an HTML report from variant data in a CSV file
def generate_report(input_csv, output_html, sample_name = None):
    df = pd.read_csv(input_csv)

    #Choosing the columns for HTML file
    base_cols = ["chrom", "pos", "rsid", "clinical_significance", "genes"]
    sample_cols = [col for col in df.columns if col not in base_cols and col not in ["ref", "alt", "condition", "disease_names"]]
    if sample_cols:
        base_cols.append(sample_cols[0])

    #Reduce DataFrame to selected columns for a clean report
    df = df[base_cols]
    table_html = df.to_html(index=False, escape=False)
    full_html = HTML_TEMPLATE.format(table=table_html)

    #Saving the final HTML report
    with open(output_html, "w") as f:
        f.write(full_html)


