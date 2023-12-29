function S_tar_cc = ClutterCancellation_Doppler(S_tar,S_ref)
S_tar=S_tar.';
S_ref=S_ref.';
S_tar_cc = S_tar-S_ref*(S_ref'*S_tar/(norm(S_ref)^2));
S_tar_cc=S_tar_cc.';
end