from cfmtoolbox import app, CFM
from cfmtoolbox_ilp_encoder.mulitsetILP import create_ilp_multiset_encoding, create_const_name


@app.command()
def encode_to_ilp_multiset(cfm: CFM) -> str:
    encoding = ""
    solver = create_ilp_multiset_encoding(cfm)
    print("Encoding")
    return solver.ExportModelAsLpFormat(False)

@app.command()
def run_ilp_solver_maximize_cardinalities(cfm: CFM):
    solver = create_ilp_multiset_encoding(cfm)
    for feature in cfm.features:
        solver.Maximize(solver.LookupVariable(create_const_name(feature)))
        solver.Solve()
        #solver.EnableOutput()
        print(print(f"{feature.name} = {solver.Objective().Value():0.1f}"))
    print(solver.ExportModelAsLpFormat(False))