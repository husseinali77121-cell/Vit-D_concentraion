import streamlit as st
import numpy as np

# إعدادات واجهة التطبيق
st.set_page_config(page_title="حاسبة التركيز الذكية", layout="centered")
st.title("🧪 نظام حساب منحنى المعايرة")
st.write("حساب التركيز باستخدام Cubic Spline")

# بيانات المعايرة الافتراضية
x_data = np.array([0.0000, 0.0104, 0.0400, 0.0798, 0.1565])
y_data = np.array([0, 5, 18.1, 41.4, 97.3])

def natural_cubic_spline(x, y):
    n = len(x) - 1
    h = np.diff(x)
    alpha = np.zeros(n)
    for i in range(1, n):
        alpha[i] = (3/h[i])*(y[i+1]-y[i]) - (3/h[i-1])*(y[i]-y[i-1])
    l, mu, z = np.ones(n+1), np.zeros(n+1), np.zeros(n+1)
    for i in range(1, n):
        l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i]
    a = y
    b, c, d = np.zeros(n), np.zeros(n+1), np.zeros(n)
    for j in range(n-1, -1, -1):
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = (y[j+1]-y[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    return a[:-1], b, c[:-1], d, x

a, b, c, d, knots = natural_cubic_spline(x_data, y_data)

# مدخلات المستخدم
abs_input = st.number_input("أدخل قيمة الامتصاصية (Absorbance):", format="%.4f")

if st.button("حساب التركيز"):
    # تحذير لو الامتصاصية برا مدى الستاندردات (الـ cubic برا المدى بيدّي أرقام غير موثوقة)
    if abs_input < knots[0] or abs_input > knots[-1]:
        st.warning(f"⚠️ الامتصاصية {abs_input:.4f} برا مدى المعايرة "
                   f"({knots[0]:.4f} – {knots[-1]:.4f}). النتيجة تقديرية وغير موثوقة؛ "
                   "يُنصح بالتخفيف وإعادة القياس.")
    idx = np.searchsorted(knots, abs_input) - 1
    idx = max(0, min(idx, len(knots)-2))
    dx = abs_input - knots[idx]
    result = a[idx] + b[idx]*dx + c[idx]*dx**2 + d[idx]*dx**3
    
    st.success(f"التركيز المحسوب هو: {result:.4f}")
    
    # عرض معلومات القطعة المستخدمة
    with st.expander("تفاصيل المعادلات"):
        st.write(f"القطعة المستخدمة: من {knots[idx]} إلى {knots[idx+1]}")
        st.code(f"A={a[idx]:.4f}, B={b[idx]:.4f}, C={c[idx]:.4f}, D={d[idx]:.4f}")
