PYTHON ?= python3

.PHONY: setup pipeline dashboard

# Install all project dependencies.
setup:
	$(PYTHON) -m pip install -r requirements.txt

# Run the full pipeline end to end with no manual intervention:
# initialize + load the database (Part 1), then generate all output
# tables and plots (Parts 2-4).
pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py

# Start the local interactive dashboard server (default http://0.0.0.0:8050).
dashboard:
	$(PYTHON) app.py
