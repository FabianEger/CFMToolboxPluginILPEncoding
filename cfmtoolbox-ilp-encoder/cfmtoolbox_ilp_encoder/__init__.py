from cfmtoolbox import app, CFM, Feature
from cfmtoolbox_ilp_encoder.mulitsetILP import create_ilp_multiset_encoding, create_const_name, \
    get_max_interval_value


@app.command()
def encode_to_ilp_multiset(cfm: CFM) -> str:
    encoding = ""
    solver = create_ilp_multiset_encoding(cfm)
    return solver.ExportModelAsLpFormat(False)

@app.command()
def run_ilp_solver_maximize_cardinalities(cfm: CFM):
    solver = create_ilp_multiset_encoding(cfm)
    find_actual_max(solver,cfm.root,1)

        #print(print(f"{feature.name} = {solver.Objective().Value():0.1f}"))
    '''
    for variable in solver.variables():
        print(f"{variable.name()} = {variable.solution_value()}")
     '''
    #print(solver.ExportModelAsLpFormat(False))

def find_actual_max(solver, feature: Feature, max_parent_cardinality: int):

    solver.Maximize(solver.LookupVariable(create_const_name(feature)))
    solver.Solve()
    #solver.EnableOutput()

    if int(solver.Objective().Value()) > 1:
        actual_max = int(solver.Objective().Value()) / max_parent_cardinality
        new_max = round(max_parent_cardinality * actual_max, None)
    else:
        actual_max = int(solver.Objective().Value())
        new_max = max_parent_cardinality
    if actual_max < get_max_interval_value(feature.instance_cardinality.intervals):
        print(feature.name + ": ")
        print("Given feature instance cardinality: " + str(get_max_interval_value(
            feature.instance_cardinality.intervals)) + "\n")
        print("Actual Max Feature Instance Cardinality " + str(round(actual_max, None)) + "\n")

    for child in feature.children:
        find_actual_max(solver, child, new_max)


    #        for variable in solver.variables():
       #     print(f"{variable.name()} = {variable.solution_value()}")


@app.command()
def run_ilp_solver_with_multisetencoding_gap_detection(cfm: CFM):
    """
    :param cfm: The input Cardinality-based Feature Model, which gets encoded
    :return: Prints the result of the SMT solver for each constant asserted to their possible cardinalities.
    """
    list_features = cfm.features
    solver = create_ilp_multiset_encoding(cfm)

    print("Searching for Gaps...")
    for feature in list_features:
        for interval in feature.instance_cardinality.intervals:
            for cardinality in range(interval.lower, interval.upper + 1):
                new_solver = create_ilp_multiset_encoding(cfm)
                constraint = new_solver.Constraint(cardinality, cardinality)
                constraint.SetCoefficient(solver.LookupVariable(create_const_name(feature)),
                                  1)
                new_solver.Maximize(solver.LookupVariable(create_const_name(feature)))
                status = new_solver.Solve()
                if status == new_solver.INFEASIBLE:
                    gap = "Gap at: " + str(cardinality) + " in Feature: " + feature.name + " "
                    print(gap)