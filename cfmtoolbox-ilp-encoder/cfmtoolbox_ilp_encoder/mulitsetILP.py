from pyexpat import features

from cfmtoolbox import CFM, Feature, Interval
from cfmtoolbox.plugins.big_m import get_global_upper_bound
from ortools.init.python import init
from ortools.linear_solver import pywraplp
from ortools.linear_solver.pywraplp import Solver
from cfmtoolbox.models import Constraint

big_M = 0

def create_ilp_multiset_encoding(cfm: CFM, maximization: bool):
    print("Encoding CFM...")

    global big_M
    big_M = get_global_upper_bound(cfm.root) * 2
    print("Big-M: " + str(get_global_upper_bound(cfm.root)))
    # Create the linear solver with the GLOP backend will give double values as result,
    # which leeds to wrong results. Therefore CBC MILP is needed
    solver = pywraplp.Solver.CreateSolver("CBC")
    if not solver:
        print("Could not create solver GLOP")
        return None

    create_ilp_multiset_variables(cfm, solver)

    if maximization:
        create_ilp_constraints_for_group_type_cardinalities(cfm.root,solver)


    create_ilp_constraints_for_group_type_cardinalities_with_less_max_than_features(cfm.root,
                                                                                    solver)

    create_ilp_constraints_for_feature_instance_cardinalities(cfm.root,solver)


    create_ilp_constraints_for_group_instance_cardinalities(cfm.root,solver)


    create_ilp_constraints(cfm.constraints,solver)

    print("Encoding complete.")

    #print("Number of variables =", solver.NumVariables())
    #print("Number of constraints =", solver.NumConstraints())
    return solver


def create_ilp_constraints_for_group_type_cardinalities(feature: Feature, solver:pywraplp.Solver):
    global big_M

    if len(feature.children) != 0:
        max_upperbound = get_max_interval_value(feature.group_type_cardinality.intervals)
        min_lowerbound = get_min_interval_value(feature.group_type_cardinality.intervals)

        constraint_lower = solver.Constraint(0, solver.infinity())
        constraint_upper = solver.Constraint(-solver.infinity(), 0)

        #constraint = solver.Constraint(min_lowerbound, max_upperbound)

        for child in feature.children:
            helper_name_1 = "Active_helper_" + child.name + "_1"
            helper = solver.BoolVar(helper_name_1)
            solver.Add(solver.LookupVariable(create_const_name(child)) - solver.LookupVariable(
                create_const_name(feature)) + big_M * (1 - helper) >= 0)

            solver.Add(solver.LookupVariable(create_const_name(child)) - solver.LookupVariable(
                create_const_name(feature)) - 1 <= big_M * helper)

            solver.Add(solver.LookupVariable(create_const_name_activ(child)) <= helper)
            solver.Add(solver.LookupVariable(create_const_name_activ(feature)) <=
                       solver.LookupVariable(create_const_name_activ_global(feature)))

            solver.Add(helper + solver.LookupVariable(create_const_name_activ_global(feature)) <=
                       1 + 2 * solver.LookupVariable(create_const_name_activ(child)))

            constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name_activ(child)), 1)
            constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name_activ(child)), 1)


        constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(feature)),
                                        -min_lowerbound)
        constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(
            feature)),
                                        -max_upperbound)

        for child in feature.children:
            create_ilp_constraints_for_group_type_cardinalities(child, solver)


def create_ilp_constraints_for_group_type_cardinalities_with_less_max_than_features(feature: Feature, solver:Solver):
    global big_M


    if len(feature.children) != 0:
        max_upperbound = get_max_interval_value(feature.group_type_cardinality.intervals)
        min_lowerbound = get_min_interval_value(feature.group_type_cardinality.intervals)

        constraint_lower = solver.Constraint(0, solver.infinity())
        constraint_upper = solver.Constraint(-solver.infinity(), 0)

        #constraint = solver.Constraint(min_lowerbound, max_upperbound)


        for child in feature.children:

            # child >= 1 - M  +  M * child_active  If child >= 1, child_active must be 1
            #helper_name_1 = "Active_helper_" + child.name + "_1"
            #solver.BoolVar(helper_name_1)

            constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(
                child)), 1)
            constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(
                child)), 1)


        constraint_lower.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(feature)),
                                        -min_lowerbound)
        constraint_upper.SetCoefficient(solver.LookupVariable(create_const_name_activ_global(
            feature)),
                                        -max_upperbound)
        #print(min_lowerbound)

        for child in feature.children:
            create_ilp_constraints_for_group_type_cardinalities_with_less_max_than_features(child, solver)


def create_ilp_constraints_for_feature_instance_cardinalities(feature_instance: Feature,
                                                              solver:Solver):
    """
        No compound intervals can be supported in the ILP encoding, which is why max and min need to be calculated
        This needs to be changed by adding new variables like in the constraint definition.
    """


    max_upperbound = get_max_interval_value(feature_instance.instance_cardinality.intervals)
    min_lowerbound = get_min_interval_value(feature_instance.instance_cardinality.intervals)


    feature_instance_variable = solver.LookupVariable(create_const_name(feature_instance))

    if len(feature_instance.instance_cardinality.intervals) > 1:

        interval_literals = []

        for i, interval in enumerate(feature_instance.instance_cardinality.intervals):
            helper_name = "helper_instance_interval_" + feature_instance.name + "_" + str(i)
            interval_helper = solver.BoolVar(helper_name)



            if feature_instance.parent is not None:
                parent_variable = solver.LookupVariable(create_const_name(feature_instance.parent))
                solver.Add(feature_instance_variable >= interval.lower * parent_variable - big_M * (1 -
                                                                                                interval_helper))

                solver.Add(feature_instance_variable <= interval.upper * parent_variable +
                                                                                big_M * (1 - interval_helper))


            else:
                solver.Add(
                    feature_instance_variable >= interval.lower - big_M * (1 -
                                                                            interval_helper))

                solver.Add(
                    feature_instance_variable <= interval.upper + big_M * (1-
                                                                          interval_helper))

            interval_literals.append(interval_helper)

        solver.Add(sum(interval_literals) == 1)

    else:
        if feature_instance.parent is not None:
            parent_variable = solver.LookupVariable(create_const_name(feature_instance.parent))

            solver.Add(feature_instance_variable >= min_lowerbound * parent_variable)
            solver.Add(feature_instance_variable <= max_upperbound * parent_variable)
        else:
            solver.Add(feature_instance_variable <= min_lowerbound)
            solver.Add(feature_instance_variable <= max_upperbound)


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
       # print(max_upperbound)
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
        create_constraint_for_intervals(solver, i, constraint.first_feature,
                                        constraint.first_cardinality.intervals,
                                        Interval(1, 3))
        create_constraint_for_intervals(solver, i, constraint.second_feature,
                                        constraint.second_cardinality.intervals,
                                        Interval(4, 6))
        helper1 = solver.LookupVariable("helper_constraint_" + str(i) +
                                                       "_" + str(2))
        helper2 = solver.LookupVariable("helper_constraint_" + str(i) + "_"
             + str(
                5))
        if not constraint.require:

            solver.Add(helper1 + helper2 <= 1)

        else:

            solver.Add(helper1 <= helper2)



def create_constraint_for_intervals(solver:Solver, constraint_number:int, feature:Feature ,
                                    cardinality: list[Interval],
                                    constants_interval: Interval):
    for i in range(constants_interval.lower,constants_interval.upper+1):
        solver.BoolVar("helper_constraint_" + str(constraint_number) + "_" + str(i))




    lower_cardinality = cardinality.__getitem__(0).lower
    upper_cardinality =  cardinality.__getitem__(0).upper
    lower_bound = lower_cardinality - 1 if lower_cardinality > 0 else 0

    helper1 = solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.lower))
    helper2 = solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.lower + 1))
    helper3 = solver.LookupVariable("helper_constraint_" + str(constraint_number) +
                                                       "_" + str(constants_interval.upper))

    feature_variable = solver.LookupVariable(create_const_name(feature))

    solver.Add(feature_variable <= lower_bound * helper1 + upper_cardinality * helper2 + big_M *
               helper3)

    solver.Add(feature_variable >=  lower_cardinality * helper2 + (upper_cardinality + 1) * helper3)


    if lower_cardinality == 0:
        excludes = solver.Constraint(1, 1)
        excludes.SetCoefficient(helper2, 1)

        excludes.SetCoefficient(helper3, 1)
    else:
        excludes = solver.Constraint(1, 1)

        excludes.SetCoefficient(helper1, 1)
        excludes.SetCoefficient(helper2, 1)

        excludes.SetCoefficient(helper3, 1)


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

        if feature.parent is None:
            max_cardinality = 1
        else:
            max_cardinality = get_max_interval_value(feature.parent.instance_cardinality.intervals)

        solver.IntVar(0, max_cardinality * big_M, create_const_name(feature)) # Big M is needed here
        # because
        # the
        # solver needs the variables to have a maximum
        solver.BoolVar(create_const_name_activ(feature))

        solver.BoolVar(create_const_name_activ_global(feature))

        # if a feature is local active it is also global active
        solver.Add(solver.LookupVariable(create_const_name_activ_global(feature)) -
                                solver.LookupVariable(create_const_name_activ(feature)) >= 0)
        eps = 1




        if feature.parent is None:
            solver.Add(solver.LookupVariable(create_const_name_activ(feature)) +
                       solver.LookupVariable(create_const_name_activ_global(feature)) == 2)
        else:
            # if a feature has >= 1 instances it is considered global activ
            solver.Add(solver.LookupVariable(
                create_const_name(feature)) <= 1 - 0.0001 + solver.LookupVariable(
                create_const_name_activ_global(
                    feature)) * big_M)
            solver.Add(solver.LookupVariable(
                create_const_name(feature)) >= 1 - (1 - solver.LookupVariable(
                create_const_name_activ_global(
                    feature))) * big_M)










def create_const_name(feature: Feature) -> str:
    return "Feature_" + feature.name

def create_const_name_activ(feature: Feature) -> str:
    return create_const_name(feature) + "_activ"

def create_const_name_activ_global(feature: Feature) -> str:
    return create_const_name(feature) + "_activ_global"