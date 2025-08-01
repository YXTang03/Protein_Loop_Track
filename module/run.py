import argparse
import json
from .main import main

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--config', 
        help='Path to the JSON config file', 
        required=True
    )
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = json.load(f)

    main(**config)

if __name__ == '__main__':
    run()
