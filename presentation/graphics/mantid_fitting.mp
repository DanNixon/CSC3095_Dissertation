input metauml;

beginfig(1);

Activity.in_workspace("Input Workspace");
Activity.model("Fit Function");
Activity.constraints("Parameter Constraints");
Activity.ties("Parameter Ties");

Activity.fit("Fit");

Activity.fit_workspace("Fit Output");
Activity.fit_params("Fitted Parameters");
Activity.cost_function_value("Cost Function Value");


leftToRight(10)(
    in_workspace,
    model);
leftToRight(10)(
    constraints,
    ties);
Group.inputs_a(
    in_workspace,
    model);
Group.inputs_b(
    constraints,
    ties);

leftToRight(10)(
    fit_workspace,
    fit_params);
Group.outputs_a(
    fit_workspace,
    fit_params);

topToBottom(10)(
    inputs_a,
    inputs_b);
Group.inputs(
    inputs_a,
    inputs_b);

topToBottom(10)(
    outputs_a,
    cost_function_value);
Group.outputs(
    outputs_a,
    cost_function_value);

topToBottom(25)(
    inputs,
    fit,
    outputs);


drawObjects(
    inputs,
    fit,
    outputs);


clink(transition)(in_workspace, fit);
clink(transition)(model, fit);
clink(transition)(constraints, fit);
clink(transition)(ties, fit);
clink(transition)(fit, fit_workspace);
clink(transition)(fit, fit_params);
clink(transition)(fit, cost_function_value);

endfig;

end
