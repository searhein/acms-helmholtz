function [S,R,N]=assemr(p,e,s,r,n,lump)
%s=stiffness
%coeff pcw constant on edges or linear
%R=robin matrix

  if nargin<6, lump=0; end

  ne=size(e,2); 
  np=size(p,2);
  if length(s)==1, s=s*ones(1,ne); end
  if length(r)==1, r=r*ones(1,ne); end
  if length(n)==1, n=n*ones(1,ne); end
%    if length(r)==1, r=r*ones(1,ne); elseif size(r,1)==np, r=pdeinedp(p,e,r); end
%    if length(n)==1, n=n*ones(1,ne); elseif size(n,1)==np, n=pdeinedp(p,e,n); end

  % geometric information
  d=p(:,e(2,:))-p(:,e(1,:)); %tangent direction 
  len=sqrt(sum(d.*d,1)); %length line segments

  %% assemble S: (s ds u,ds v)_Gamma 
  if length(s)==np
      s=(s(e(1,:))+s(e(2,:)))/2;
  end
  slen=s./len; 
  elmats=[ slen; -slen; ...
          -slen;  slen]; %elements matrices 1/h -1/h 
  ii=[e(1,:);e(1,:); e(2,:);e(2,:)]; %entries based on global matrix
  jj=[e(1,:);e(2,:); e(1,:);e(2,:)];
  S=sparse(ii(:),jj(:),elmats,np,np);
  
  
  %% assemble R: (cr u, v)_Gamma
  if length(r)==0
     R=sparse(np,np);
  elseif size(r,2)==ne
    if lump
      rlen=r.*len/2;
      elmats=[  rlen; 0*rlen; ...
              0*rlen;   rlen];
    else
      r1len=r.*len/6; r2len=2*r1len;
      elmats=[r2len;r1len; ...
              r1len;r2len];
    end
    ii=[e(1,:);e(1,:); e(2,:);e(2,:)];
    jj=[e(1,:);e(2,:); e(1,:);e(2,:)];
    R=sparse(ii(:),jj(:),elmats,np,np);
  elseif size(r,1)==np %% c is pw linear!
    if lump
      r1len=r(e(1,:))'.*len/2; r2len=r(e(2,:))'.*len/2;
      elmats=[  r1len; 0*r1len; ...
              0*r2len;   r2len];
    else
      r1=r(e(1,:))'; r2=r(e(2,:))';
      r1len12==r1.*len./12; r1len4=3*r1len12;
      r2len12==r2.*len./12; r2len4=3*r2len12;
      elmats=[r1len4+r2len12; r1len12+r2len12; ...
              r1len12+r2len12; r1len12+r2len4];
    end
    ii=[e(1,:);e(1,:); e(2,:),e(2,:)];
    jj=[e(1,:);e(2,:); e(1,:),e(2,:)];    
    R=sparse(ii(:),jj(:),elmats,np,np);
  else
    error('dimension mismatch')
  end

  
  % assemble N: (n, v)_Gamma
  if size(n,1)==0 || nargout<2
     N=zeros(np,1);
  elseif size(n,2)==ne
    nlen2=n.*len*0.5;
    elvecs=[nlen2;nlen2];
    i=e(1:2,:);
    N=full(sparse(i(:),ones(size(i(:))),elvecs,np,1));
  elseif size(n,1)==np
    n1=n(e(1,:))'; n2=n(e(2,:))'; 
    len6=len/6; len3=len/3; 
    elvec=[n1.*len3+n2.*len6; n1.*len6+n2.*len3];
    i=e(1:2,:);
    N=full(sparse(i(:),ones(size(i(:))),elvec,np,1));
  else
   error('size of parameter n does not match') 
  end

