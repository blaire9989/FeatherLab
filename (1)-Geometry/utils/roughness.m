function perturb = roughness(tau, xMin, xMax, yMin, yMax, N, permute)

    % Generate a 1D roughness distribution, used for simulating spatially varying barbule features
    if rem(N, 2) == 1
        N = N + 1;
    end
    frange = 1 / ((xMax - xMin) / N);
    fMin = -0.5 * frange;
    fMax = 0.5 * frange;
    freqlim = linspace(fMin, fMax, N + 1);
    amplitude = exp(-pi * pi * freqlim.^2 * tau^2);
    phase = rand(1, N / 2) * 2 * pi - pi;
    phase = [phase 0 -flip(phase)];
    freq = transpose(amplitude .* exp(1j * phase));
    [perturb, ~, ~] = ift(freq(1:N), fMin, fMax);
    perturb = real(perturb);
    if permute
        [~, index] = min(perturb);
        if index > 1
            perturb = [perturb(index : end); perturb(1 : index - 1)];
        end
    end
    rMin = min(perturb);
    rMax = max(perturb);
    perturb = perturb / (rMax - rMin) * (yMax - yMin);
    perturb = perturb - min(perturb) + yMin;

end

function [g, xMin, xMax] = ift(G, f0, fN)

    % 1D inverse Fourier transform
    N = size(G, 1);
    delta_f = (fN - f0) / N;
    Xs = 1 / delta_f;
    xMin = -0.5 * Xs;
    xMax = 0.5 * Xs;
    temp_g = ifft(G);
    N0 = N / 2;
    g = zeros(N, 1);
    inds = transpose((1:N) - N0 - 1) / N;
    g(1:N0) = (fN - f0) .* exp(2 * pi * 1j * inds(1:N0) * Xs * f0) .* temp_g((N0 + 1):N);
    g((N0 + 1):N) = (fN - f0) .* exp(2 * pi * 1j * inds((N0 + 1):N) * Xs * f0) .* temp_g(1:N0);

end