'''Constants Used in iBuilderSB'''

import os

PROJECT_NAME = "iBuilderSB"
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
for _ in range(2):
    PROJECT_DIR = os.path.dirname(PROJECT_DIR)
    if (PROJECT_NAME in PROJECT_DIR) and (not "src" in PROJECT_DIR):
        break
DATA_DIR = os.path.join(PROJECT_DIR, "data")
REACTIONS_PATH = os.path.join(DATA_DIR, "rhea_reactions.csv")