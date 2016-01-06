function [L,R,P,x] = my_factorization_PGBC(A,b,mode)

% Input
%		A	  : nxn square matrix to be factorized as LR
%		b     : nx1 column vector such that Ax=b
%		mode  : 'WOP' - without pivoting; 'WP' - with pivoting
% Output
% 		L	  : lower triangular matrix such that LR=PA
%		R	  : upper triangular matrix such that LR=PA 
%		P     : permutation matrix such that LR=PA
%			  : P is the nxn identity matrix for mode WOP
%		x     : nx1 column vector which satisfies Ax=b

dim=length(A);
% Initialize the solution column vector.
x=zeros(dim,1);
fail_flag=false;

if strcmp(mode,'WOP')
	L=eye(dim);
	R=A;
	P=eye(dim);
	for i=1:dim-1
		% Check if the entries in column i from row i+1 downward are 0.
		% If this is the case, we can proceed to the next iteration.
		if sum(R(i+1:dim,i) == zeros(dim-i,1)) == dim-i
			continue
		end
		pivot=R(i,i);
		% Check if the pivot element of some row is 0.
		% If this is the case, set fail_flag to true.
		if pivot==0
			fprintf('Mode WOP failed: pivot for row %d is 0. \n', i)
			fprintf('Switched to mode WP.\n')
			fail_flag=true;
			break
		end
		L_i=eye(dim);
		L_i_inv=eye(dim);
		% Construct the ith Frobenius matrix and its inverse.
		for j=i+1:dim
			L_i(j,i)=-R(j,i)/pivot;
			L_i_inv(j,i)=-L_i(j,i);
		end
		% Construct the triangular matrices L and R.
		L=L*L_i_inv;
		R=L_i*R;
	end
	if sum(sum(L*R == A)) ~= dim*dim
		fprint('Mode WOP failed: LR not equal to A. \n')
	    fprintf('Switched to mode WP.\n')
		fail_flag = true;
	end
	if fail_flag == false
		fprintf('Mode WOP is successful.\n');
	end
end

if strcmp(mode,'WP') || fail_flag==true
	L=eye(dim);
	R=A;
	P=eye(dim);
	I=eye(dim);
	for i=1:dim-1
		% Obtain the index of the maximum element among the
		% elements in column i and rows i to dim.
		[~,rel_max_index]=max(abs(R(i:dim,i)));
		% Compute the row number of the maximum element.
		max_index=rel_max_index+i-1;
		P_i=eye(dim);
		% Construct the ith permutation matrix.
		P_i(i,1:dim)=I(max_index,1:dim);
		P_i(max_index,1:dim)=I(i,1:dim);
		% Bring the maximum element to the top.
		R=P_i*R;
		pivot=R(i,i);
		L_i=eye(dim);
		L_i_inv=eye(dim);
		% Construct the ith Frobenius matrix and its inverse.
		for j=i+1:dim
			L_i(j,i)=-R(j,i)/pivot;
			L_i_inv(j,i)=-L_i(j,i);
		end;
		% Construct the permutation matrix P and the 
		% triangular matrices L and R.
		P=P_i*P;
		R=L_i*R;
		% For any permutation matrix P, P_inv = P'
		L=L*P_i'*L_i_inv;
	end		
	% Premultiply L by P to cancel out P_inv in L.
	L=P*L;
	fprintf('Mode WP is successful.\n');
% Check if mode entered is neither WOP nor WP.
elseif ~strcmp(mode,'WOP')
	fprintf('Invalid mode. WOP : w/o pivoting, WP : w/ pivoting.\n')
	return
end

% Perform forward substitution on the system Lz=Pb to solve for z.
b=P*b;
z=zeros(dim,1);
z(1)=b(1)/L(1,1);
for i = 2:dim
	z(i)=(b(i)-L(i,1:i-1)*z(1:i-1))/L(i,i);
end

% Perform backward substitution on the system Rx=z to solve for x.
x(dim)=z(dim)/R(dim,dim);
for i = dim-1:-1:1
	x(i)=(z(i)-R(i,i+1:dim)*x(i+1:dim))/R(i,i);
end

%Prepared by Pio Gabrielle B Calderon 2010-09736