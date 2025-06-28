# OAS Database Processing

This repo contains code for downloading and post-processing the OPIG OAS database.


# Setup

```
python -m venv venv
source venv/bin/activate
python -m pip install -r requirements.txt
```

# Usage

Downloads the OAS database csv files in `oas_paired_urls.txt` and extracts the metadata header line.

```bash
./download-oas && ./extract-json-headers
```

Converts the downloaded csv files to parquet files.
```bash
./csv2parquet.py $HOME/Downloads/oas_paired outputdir commonheaders
```
