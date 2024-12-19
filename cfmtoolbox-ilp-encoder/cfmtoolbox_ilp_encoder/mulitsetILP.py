from cfmtoolbox import CFM, Feature, Interval
from ortools.init.python import init
from ortools.linear_solver import pywraplp
from ortools.linear_solver.pywraplp import Solver, Constraint, ExportModelAsLpFormat


def create_ilp_multiset_encoding(cfm: CFM):
    print("Encoding CFM...")

    # Create the linear solver with the GLOP backend.
    solver = pywraplp.Solver.CreateSolver("GLOP")
    if not solver:
        print("Could not create solver GLOP")
        return None

    create_ilp_multiset_variables(cfm, solver)
    print("Number of variables =", solver.NumVariables())

    create_ilp_constraints_for_group_type_cardinalities(cfm.root,solver)
    print("Number of constraints =", solver.NumConstraints())

    create_ilp_constraints_for_feature_instance_cardinalities(cfm.root,solver)
    print("Number of constraints =", solver.NumConstraints())

    create_ilp_constraints_for_group_instance_cardinalities(cfm.root,solver)
    print("Number of constraints =", solver.NumConstraints())

    return solver

def create_ilp_constraints_for_group_type_cardinalities(feature: Feature, solver:Solver):

    if len(feature.children) != 0:
        max_upperbound = get_max_interval_value(feature.group_type_cardinality.intervals)
        min_lowerbound = get_min_interval_value(feature.group_type_cardinality.intervals)

        constraint_lower = solver.Constraint(0, solver.infinity())
        constraint_upper = solver.Constraint(-solver.infinity(), 0)

        #constraint = solver.Constraint(min_lowerbound, max_upperbound)

        for child in feature.children:
            constraint_lower.SetCoefficient(solver.LookupVariable(creat_const_name_activ(child)),1)
            constraint_upper.SetCoefficient(solver.LookupVariable(creat_const_name_activ(child)),1)

        constraint_lower.SetCoefficient(solver.LookupVariable(creat_const_name_activ(feature)),
                                        -min_lowerbound)
        constraint_upper.SetCoefficient(solver.LookupVariable(creat_const_name_activ(feature)),
                                        -max_upperbound)

        for child in feature.children:
            create_ilp_constraints_for_group_type_cardinalities(child, solver)


def create_ilp_constraints_for_feature_instance_cardinalities(feature_instance: Feature,
                                                              solver:Solver):
    """
        No compound intervals can be supported in the ILP encoding, which is why max and min need to be calculated
    """
    max_upperbound = get_max_interval_value(feature_instance.instance_cardinality.intervals)
    min_lowerbound = get_min_interval_value(feature_instance.instance_cardinality.intervals)



    if feature_instance.parent is not None:
        constraint_lower = solver.Constraint(0, solver.infinity())
        constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)),1)
        constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance.parent)),-min_lowerbound)

        constraint_upper = solver.Constraint(-solver.infinity(), 0)
        constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)), 1)
        constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance.parent)), -max_upperbound)
    else:
        constraint = solver.Constraint(min_lowerbound, max_upperbound)
        constraint.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)),1)



    for child in feature_instance.children:
        create_ilp_constraints_for_feature_instance_cardinalities(child, solver)



def create_ilp_constraints_for_group_instance_cardinalities(feature_instance: Feature,
                                                              solver:Solver):
    """
        No compound intervals can be supported in the ILP encoding, which is why max and min need to be calculated
    """
    if len(feature_instance.children) != 0:
        max_upperbound = get_max_interval_value(feature_instance.group_instance_cardinality.intervals)
        min_lowerbound = get_min_interval_value(feature_instance.group_instance_cardinality.intervals)

        constraint_lower = solver.Constraint(0, solver.infinity())
        constraint_upper = solver.Constraint(-solver.infinity(), 0)

        for child in feature_instance.children:
            constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name(child)),1)
            constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name(child)),1)

        constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)),-min_lowerbound)

        constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)),-max_upperbound)

        for child in feature_instance.children:
            create_ilp_constraints_for_group_instance_cardinalities(child, solver)






def get_max_interval_value(intervals: list[Interval])-> int:
    if len(intervals) == 0:
        return 0
    else:
        max_upperbound = intervals[0].upper
        for interval in intervals:
            if interval.upper > max_upperbound:
                max_upperbound = interval.upper
        return max_upperbound


def get_min_interval_value(intervals: list[Interval])-> int:
    if len(intervals) == 0:
        return 0
    else:
        min_lowerbound = intervals[0].lower
        for interval in intervals:
            if interval.lower < min_lowerbound:
                min_lowerbound = interval.lower
        return min_lowerbound


def create_ilp_multiset_variables(cfm: CFM, solver: Solver):
    for feature in cfm.features:
        solver.NumVar(0,  10000000, create_const_name(feature)) # Big M is needed here
        solver.NumVar(0, 1, creat_const_name_activ(feature))

        constraint = solver.Constraint(-solver.infinity(), 0)
        constraint.SetCoefficient(solver.LookupVariable(creat_const_name_activ(feature)),1)
        constraint.SetCoefficient(solver.LookupVariable(create_const_name(feature)), -1)

def create_const_name(feature: Feature) -> str:
    return "Feature_" + feature.name

def creat_const_name_activ(feature: Feature) -> str:
    return create_const_name(feature) + "_activ"