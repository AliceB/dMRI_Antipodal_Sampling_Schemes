function [ q_t theta_t ] = quad_wts_ssht( L_h )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% same as jason's paper, here L_h is L-1 or L = L_h+1


theta_t =pi/(2*L_h+1):2*pi/(2*L_h+1):pi;

v_t1 = 0;
for mm=-L_h:1:L_h
v_t1 = v_t1 + reflected_wts(-mm)*exp(1i*mm*theta_t);


end


v_t2 = 0;
for mm=-L_h:1:L_h
v_t2 = v_t2 + reflected_wts(-mm)*exp(-1i*mm*theta_t);


end


q_t = real(4*pi*pi*(v_t1+v_t2)/(2*L_h+1));

% for pi
q_t(L_h+1) = real(4*pi*pi*v_t1(L_h+1)/(2*L_h+1));


end