function q = interface_flux(~,CB,CE,params)

    % Same interface flux used by both brine and EPS.
    % Positive q means nutrient leaves brine and enters EPS.

    gamma = params.gamma;

    % simple closure for the shared interface flux
    q = gamma*(CB - CE);

end