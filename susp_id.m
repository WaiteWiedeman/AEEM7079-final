function dx=susp_id(x,k_t,k_s,m_b,m_w,wk)
dx = zeros(5,1);
dx(1) = x(3) - x(4);
dx(2) = x(4) - wk;
dx(3) = -k_s/m_b*x(1) - x(5)/m_b*x(3) + x(5)/m_b*x(4);
dx(4) = -k_s/m_w*x(1) - k_t/m_w*x(2) + x(5)/m_w*x(3) - x(5)/m_w*x(4);
dx(5) = 0;