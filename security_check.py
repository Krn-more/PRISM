import os
import sys
import datetime
import requests
import tkinter as tk
from tkinter import messagebox

# Configuration
AUTHORIZED_DOMAIN = "CODE1"
INTERNAL_URL = "https://share.philips.com"
EXPIRATION_DATE = datetime.date(2026, 9, 30)

def show_error_and_exit(title, message):
    """Shows a native Windows error dialog and exits the application."""
    root = tk.Tk()
    root.withdraw() # Hide the main tk window
    
    # Try to keep it on top
    root.attributes('-topmost', True)
    
    messagebox.showerror(title, message)
    sys.exit(1)

def check_domain():
    """Check if the user is on the authorized active directory domain."""
    user_domain = os.environ.get("USERDOMAIN", "")

    if user_domain.upper() != AUTHORIZED_DOMAIN.upper():
        show_error_and_exit(
            "Security Verification Failed",
            f"This application is restricted to the {AUTHORIZED_DOMAIN} network domain.\n\n"
            f"Current Domain: {user_domain if user_domain else 'Unknown'}\n"
            "Access Denied."
        )

def check_network():
    """Check if we can ping/reach the internal company network."""
    try:
        # A lightweight HEAD request with a short timeout
        response = requests.head(INTERNAL_URL, timeout=3.0)
        # We don't necessarily care if it's a 200 or 401, as long as it resolved and responded
    except requests.RequestException:
        show_error_and_exit(
            "Network Verification Failed",
            f"Unable to reach the internal company network ({INTERNAL_URL}).\n\n"
            "Please ensure you are connected to the company network or VPN and try again."
        )

def check_expiration():
    """Check if the hardcoded time bomb has expired."""
    current_date = datetime.date.today()
    if current_date > EXPIRATION_DATE:
        show_error_and_exit(
            "License Expired",
            f"This version of the application expired on {EXPIRATION_DATE.strftime('%B %d, %Y')}.\n\n"
            "Please contact the IT department or application maintainer for an updated version."
        )

def run_security_checks():
    """Run all security checks sequentially."""
    check_expiration()
    check_domain()
    check_network()

if __name__ == "__main__":
    # For testing the script directly
    run_security_checks()
    print("All security checks passed.")
