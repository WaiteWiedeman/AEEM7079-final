function dx=susp(x,c,k_t,k_s,m_b,m_w,wk)
dx = zeros(4,1);
dx(1) = x(3) - x(4);
dx(2) = x(4) - wk;
dx(3) = -k_s/m_b*x(1) - c/m_b*x(3) + c/m_b*x(4);
dx(4) = -k_s/m_w*x(1) - k_t/m_w*x(2) + c/m_w*x(3) - c/m_w*x(4);
