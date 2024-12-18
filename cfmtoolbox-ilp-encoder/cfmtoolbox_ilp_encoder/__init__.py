from cfmtoolbox import app, CFM
from cfmtoolbox_ilp_encoder.mulitsetILP import create_ilp_multiset_encoding


@app.command()
def encode_to_ilp_multiset(cfm: CFM) -> str:
    encoding = ""
    create_ilp_multiset_encoding(cfm)
    print("Encoding")
    return encoding