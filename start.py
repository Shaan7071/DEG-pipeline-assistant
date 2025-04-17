import os
import subprocess

# Get the PORT environment variable, defaulting to 8501 if not set
port = os.environ.get('PORT', '8501')

# Start Streamlit with the correct port
subprocess.run(['streamlit', 'run', 'app_v2.py', 
                '--server.port=' + port, 
                '--server.address=0.0.0.0'])
