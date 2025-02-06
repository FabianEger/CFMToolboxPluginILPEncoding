from pyexpat import features

from cfmtoolbox import CFM, Feature, Interval
from cfmtoolbox.plugins.big_m import get_global_upper_bound
from ortools.init.python import init
from ortools.linear_solver import pywraplp
from ortools.linear_solver.pywraplp import Solver
from cfmtoolbox.models import Constraint

big_M = 0

def create_ilp_multiset_encoding(cfm: CFM):
    print("Encoding CFM...")

    global big_M
    big_M = get_global_upper_bound(cfm.root)

    # Create the linear solver with the GLOP backend will give double values as result,
    # which leeds to wrong results. Therefore CBC MILP is needed
    solver = pywraplp.Solver.CreateSolver("CBC")
    if not solver:
        print("Could not create solver GLOP")
        return None

    create_ilp_multiset_variables(cfm, solver)


    create_ilp_constraints_for_group_type_cardinalities(cfm.root,solver)


    create_ilp_constraints_for_feature_instance_cardinalities(cfm.root,solver)


    create_ilp_constraints_for_group_instance_cardinalities(cfm.root,solver)


    create_ilp_constraints(cfm.constraints,solver)

    print("Encoding complete.")

    #print("Number of variables =", solver.NumVariables())
    #print("Number of constraints =", solver.NumConstraints())
    return solver


def create_ilp_constraints_for_group_type_cardinalities(feature: Feature, solver:Solver):
    global big_M
    if feature.parent is None:
        constraint = solver.Constraint(1, 1)
        constraint.SetCoefficient(solver.LookupVariable(create_const_name_activ(feature)),
                                  1)

    if len(feature.children) != 0:
        max_upperbound = get_max_interval_value(feature.group_type_cardinality.intervals)
        min_lowerbound = get_min_interval_value(feature.group_type_cardinality.intervals)

        constraint_lower = solver.Constraint(0, solver.infinity())
        constraint_upper = solver.Constraint(-solver.infinity(), 0)

        #constraint = solver.Constraint(min_lowerbound, max_upperbound)


        for child in feature.children:


            # child >= parent - M  +  M * child_active  If child >= parent, child_active must be 1
            #helper_name_1 = "Active_helper_" + child.name + "_1"
            #solver.BoolVar(helper_name_1)

            constraint_upper_local = solver.Constraint(-solver.infinity(), big_M)
            constraint_upper_local.SetCoefficient(solver.LookupVariable(create_const_name(
                feature)), 1)
            constraint_upper_local.SetCoefficient(solver.LookupVariable(create_const_name(
                child)), -1)
            constraint_upper_local.SetCoefficient(solver.LookupVariable(create_const_name_activ(
                child)), big_M)

            # child <= parent - 1 + M * helper1  If child < parent, helper1 must be 0
            constraint_lower_local = solver.Constraint(1, solver.infinity())
            constraint_lower_local.SetCoefficient(solver.LookupVariable(create_const_name(
                feature)),1)
            constraint_lower_local.SetCoefficient(solver.LookupVariable(
                create_const_name(
                child)),-1)
            constraint_lower_local.SetCoefficient(solver.LookupVariable(create_const_name_activ(
                child)),big_M)


            constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name_activ(child)), 1)
            constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name_activ(child)), 1)


        constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(feature)),
                                        -min_lowerbound)
        constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(
            feature)),
                                        -max_upperbound)

        for child in feature.children:
            create_ilp_constraints_for_group_type_cardinalities(child, solver)


def create_ilp_constraints_for_feature_instance_cardinalities(feature_instance: Feature,
                                                              solver:Solver):
    """
        No compound intervals can be supported in the ILP encoding, which is why max and min need to be calculated
        This needs to be changed by adding new variables like in the constraint definition.
    """

    ''' 
    max_upperbound = get_max_interval_value(feature_instance.instance_cardinality.intervals)
    min_lowerbound = get_min_interval_value(feature_instance.instance_cardinality.intervals)

    if feature_instance.parent is not None:
        constraint_lower = solver.Constraint(0, solver.infinity())
        constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)), 1)
        constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance.parent)), -min_lowerbound)

        constraint_upper = solver.Constraint(-solver.infinity(), 0)
        constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)), 1)
        constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance.parent)), -max_upperbound)
    else:
        constraint = solver.Constraint(min_lowerbound, max_upperbound)
        constraint.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)), 1)

    for child in feature_instance.children:
        create_ilp_constraints_for_feature_instance_cardinalities(child, solver)




    '''

    only_one_interval_constraint = solver.Constraint(0, 0)


    if feature_instance.parent is not None:
        only_one_interval_constraint.SetCoefficient(solver.LookupVariable(
            create_const_name_activ_global(feature_instance.parent)), -1)
        constraint_lower = solver.Constraint(-solver.infinity(),0)
        constraint_upper = solver.Constraint(0,solver.infinity())
        max_upperbound = get_max_interval_value(feature_instance.parent.instance_cardinality.intervals)
        min_lowerbound = get_min_interval_value(
            feature_instance.parent.instance_cardinality.intervals)
    else:
        only_one_interval_constraint.SetCoefficient(solver.LookupVariable(
            create_const_name_activ_global(feature_instance)), -1)
        constraint = solver.Constraint(0, 0)

    # Feature >= parent * interval.lower
    #  Feature >= (helper * firstinterval.lower * parent) + (helper2 * secondinterval.lower *
    #  parent)
    # Feature <= helper * firstinterval + helper2 * secondinteval




    for i,interval in enumerate(feature_instance.instance_cardinality.intervals):
        helper_name = "helper_instance_interval_" + feature_instance.name + "_" + str(i)
        solver.BoolVar(helper_name)
        only_one_interval_constraint.SetCoefficient(solver.LookupVariable(helper_name), 1)

        if feature_instance.parent is not None:
            # Auxiliary variable
            v = solver.IntVar(0, big_M, 'v_' +  feature_instance.name + "_" + str(i))  # v = parent
            # * helper

            # Step 1: Linearize v = parent * helper (McCormick Envelopes)
            solver.Add(v >= 0)
            solver.Add(v >= solver.LookupVariable(
                create_const_name(feature_instance.parent)) - big_M * solver.LookupVariable(
                helper_name) - big_M)
            solver.Add(v <= solver.LookupVariable(create_const_name(feature_instance.parent)))
            solver.Add(v <= big_M * solver.LookupVariable(helper_name))

            constraint_lower.SetCoefficient(v,interval.lower)
            constraint_upper.SetCoefficient(v,interval.upper)
        else:
            constraint.SetCoefficient(solver.LookupVariable(helper_name),interval.lower)

    if feature_instance.parent is not None:

        constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)),-1)
       # constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name(
        #    feature_instance.parent)),-min_lowerbound)


        constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)), -1)
        #constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name(
         #   feature_instance.parent)), -max_upperbound)
    else:
        constraint.SetCoefficient(solver.LookupVariable(create_const_name(
            feature_instance)), -1)


    for child in feature_instance.children:
        create_ilp_constraints_for_feature_instance_cardinalities(child, solver)


def create_ilp_constraints_for_group_instance_cardinalities(feature_instance: Feature,
                                                              solver:Solver):
    """
        No compound intervals are supported in the ILP encoding for group instances, which is why
        max and min
        need to be calculated
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




def create_ilp_constraints(constraints: list[Constraint], solver:Solver):
    """
    Assumption that constraints do not have compound intervals. And in the intervals are no *
    operator.
    :param constraints:
    :param solver:
    :return:
    """
    for i,constraint in enumerate(constraints):
        if not constraint.require:
            create_constraint_for_intervals(solver,i,constraint.first_feature,
                                            constraint.first_cardinality.intervals,
                                            Interval(1,3))
            create_constraint_for_intervals(solver,i,constraint.second_feature,
                                            constraint.second_cardinality.intervals,
                                            Interval(4,6))

            # the features can not both be in their intervals so the addition of the variables
            # needs to be less equal 1
            constraint = solver.Constraint(-solver.infinity(), 1)
            constraint.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(i) +
                                                       "_" + str(2)),1)
            constraint.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(i) + "_"
             + str(
                5)),1)
        else:
            # if the first feature is in the given Interval than the second needs to be in the
            # given interval -> HelperIntervalConstFeature2 >= HelperIntervalConstFeature1 -> 0
            # >= HelperIntervalConstFeature1 - HelperIntervalConstFeature2
            create_constraint_for_intervals(solver,i,constraint.first_feature,
             constraint.first_cardinality.intervals,
                                            Interval(1,3))
            create_constraint_for_intervals(solver,i,constraint.second_feature,
                                            constraint.second_cardinality.intervals,Interval(4,6))
            constraint = solver.Constraint(-solver.infinity(), 0)
            constraint.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(i) +
                                                       "_" + str(2)),1)
            constraint.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(i) +
                                                       "_" + str(5)),-1)


def create_constraint_for_intervals(solver:Solver, constraint_number:int, feature:Feature ,
                                    cardinality: list[Interval],
                                    constants_interval: Interval):
    for i in range(constants_interval.lower,constants_interval.upper+1):
        solver.BoolVar("helper_constraint_" + str(constraint_number) + "_" + str(i))

    '''
    global_active_feature_var_name = "helper_constraint_feature_global_active" + "_" + str(
        feature.name) + "_" + str(constraint_number)
    solver.BoolVar(global_active_feature_var_name)

    # x <= M * z   If z = 0, x <= 0
    #x >= epsilon - M * (1 - z)  # If z = 1, x >= epsilon
    
    global_constraint_upper = solver.Constraint(-solver.infinity(), big_M - 1)
    global_constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name(feature)), -1)
    global_constraint_upper.SetCoefficient(solver.LookupVariable(global_active_feature_var_name), big_M)
    global_constraint_lower = solver.Constraint(0, solver.infinity())
    global_constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name(feature)), -1)
    global_constraint_lower.SetCoefficient(solver.LookupVariable(global_active_feature_var_name), big_M)
    
    '''


    exclude_upper = solver.Constraint(0, solver.infinity())
    exclude_upper.SetCoefficient(solver.LookupVariable(create_const_name(
        feature)), -1)

    lower_cardinality = 0 if cardinality.__getitem__(0).lower == 0 else (cardinality.__getitem__(
        0).lower - 1)
    exclude_upper.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.lower)),
                                 lower_cardinality)
    exclude_upper.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.lower + 1)),
                                 cardinality.__getitem__(
                                     0).upper)
    exclude_upper.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.upper)), big_M)

    exclude_lower = solver.Constraint(-solver.infinity(),0)
    exclude_lower.SetCoefficient(solver.LookupVariable(create_const_name(
        feature)), -1)
    exclude_lower.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(feature)), -1)
    exclude_lower.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.lower)),1
                                 )
    exclude_lower.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.lower + 1)),
                                 cardinality.__getitem__(0).lower)
    exclude_lower.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.upper)),
                                 cardinality.__getitem__(
                                     0).upper + 1)





    excludes = solver.Constraint(0,0)
    excludes.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.lower)), 1)
    excludes.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.lower + 1)), 1)
    if not ( cardinality.__getitem__(0).lower ==  cardinality.__getitem__(0).upper):
        excludes.SetCoefficient(solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.upper)), 1)
    excludes.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(feature)),
                            -1)


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
    global big_M
    for feature in cfm.features:
        solver.IntVar(0, big_M , create_const_name(feature)) # Big M is needed here because the
        # solver needs the variables to have a maximum
        solver.IntVar(0, 1, create_const_name_activ(feature))

        solver.IntVar(0, 1, create_const_name_activ_global(feature))


        constraint = solver.Constraint(0, solver.infinity())
        constraint.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(feature)), 1)
        constraint.SetCoefficient(solver.LookupVariable(create_const_name_activ(feature)), -1)








def create_const_name(feature: Feature) -> str:
    return "Feature_" + feature.name

def create_const_name_activ(feature: Feature) -> str:
    return create_const_name(feature) + "_activ"

def create_const_name_activ_global(feature: Feature) -> str:
    return create_const_name(feature) + "_activ_global"