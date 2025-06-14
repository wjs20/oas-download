# OAS Database Processing

This repo contains code for downloading and post-processing the OPIG OAS database.


# Setup

```
python -m venv venv
source venv/bin/activate
python -m pip install -r requirements.txt
```

# Usage
```bash
./download-oas
```

```bash
./csv2parquet.py $HOME/Downloads/oas_paired outputdir commonheaders
```
