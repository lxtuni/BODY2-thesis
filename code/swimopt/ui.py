# -*- coding: utf-8 -*-
"""
swimopt 控制面板 —— 把导入模型/切换/设置搜索空间/寻优/可视化 全装进一个窗口.
运行: python ui.py    (或双击 8_控制面板.bat)
"""
import os, re, sys, json, glob, time, queue, threading, subprocess
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import tkinter.font as tkfont

# Windows 高分屏: 开启 DPI 感知, 否则字会糊
try:
    from ctypes import windll
    try: windll.shcore.SetProcessDpiAwareness(2)      # per-monitor v2
    except Exception: windll.user32.SetProcessDPIAware()
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
PY = sys.executable


class App:
    def __init__(self, root):
        self.root = root
        root.title("swimopt 控制面板 —— 四足水下机器人 最优泳姿寻优")
        root.geometry("1340x920")
        self.proc = None
        self.q = queue.Queue()
        self.hist = []            # (eval, best_fitness)
        self.cfg_path = "config.json"
        self._setup_fonts()
        self._build()
        self._refresh_models()
        self._load_cfg_into_ui()
        for e in (self.e_wyaw, self.e_wene, self.e_amax, self.e_f0, self.e_f1, self.e_omax):
            e.bind("<FocusOut>", lambda ev: self._update_dim())
            e.bind("<Return>", lambda ev: self._update_dim())
        self.root.after(300, self._pump)

    def _setup_fonts(self):
        """统一字体: 中文用微软雅黑, 数字表格用等宽, 放大字号避免糊."""
        fam = "Microsoft YaHei UI"
        avail = set(tkfont.families())
        if fam not in avail:
            fam = "Microsoft YaHei" if "Microsoft YaHei" in avail else "TkDefaultFont"
        mono = "Consolas" if "Consolas" in avail else "Courier New"
        for name, size in [("TkDefaultFont", 11), ("TkTextFont", 11), ("TkMenuFont", 11),
                           ("TkHeadingFont", 11), ("TkTooltipFont", 10)]:
            try:
                f = tkfont.nametofont(name); f.configure(family=fam, size=size)
            except Exception:
                pass
        st = ttk.Style()
        try: st.theme_use("clam")
        except Exception: pass
        st.configure(".", font=(fam, 11))
        st.configure("TLabelframe.Label", font=(fam, 11, "bold"), foreground="#1a3c6e")
        st.configure("TButton", font=(fam, 11), padding=5)
        st.configure("Treeview", font=(mono, 11), rowheight=26)
        st.configure("Treeview.Heading", font=(fam, 11, "bold"))
        self.font_mono = (mono, 11)
        self.font_big = (fam, 12, "bold")

    # ---------------- 界面 ----------------
    def _build(self):
        pad = dict(padx=6, pady=4)

        # === 顶部: 模型 ===
        f1 = ttk.LabelFrame(self.root, text=" 1. 模型 ")
        f1.pack(fill="x", **pad)
        ttk.Label(f1, text="当前模型:").pack(side="left", padx=6)
        self.cb_model = ttk.Combobox(f1, width=34, state="readonly")
        self.cb_model.pack(side="left", padx=4)
        self.cb_model.bind("<<ComboboxSelected>>", lambda e: self._switch_model())
        ttk.Button(f1, text="导入 URDF…", command=self._import_urdf).pack(side="left", padx=4)
        ttk.Button(f1, text="刷新列表", command=self._refresh_models).pack(side="left", padx=4)
        self.lb_model = ttk.Label(f1, text="", foreground="#25507a")
        self.lb_model.pack(side="left", padx=12)

        # === 搜索空间 ===
        f2 = ttk.LabelFrame(self.root, text=" 2. 搜索空间 —— 决定优化器在随机/搜索哪些参数 ")
        f2.pack(fill="x", **pad)
        self.v_freq = tk.BooleanVar(value=True)
        self.v_duty = tk.BooleanVar(value=True)
        self.v_amp = tk.BooleanVar(value=True)
        self.v_phase = tk.BooleanVar(value=True)
        self.v_off = tk.BooleanVar(value=True)
        for txt, var in [("频率 f", self.v_freq), ("占空比(发力%)", self.v_duty), ("幅度 A", self.v_amp),
                         ("相位 φ (=步态)", self.v_phase), ("偏置 offset", self.v_off)]:
            ttk.Checkbutton(f2, text=txt, variable=var,
                            command=self._update_dim).pack(side="left", padx=8)
        ttk.Label(f2, text="   相位冻结时固定为:").pack(side="left")
        self.cb_preset = ttk.Combobox(f2, width=10, state="readonly",
                                      values=["diag 对角", "fb 前后", "lr 左右", "wave 行波", "inphase 同相"])
        self.cb_preset.current(0)
        self.cb_preset.pack(side="left", padx=4)
        self.cb_preset.bind("<<ComboboxSelected>>", lambda e: self._update_dim())
        self.lb_dim = ttk.Label(f2, text="", foreground="#a0410d", font=self.font_big)
        self.lb_dim.pack(side="left", padx=14)

        # === 运行设置 ===
        f3 = ttk.LabelFrame(self.root, text=" 3. 运行设置 ")
        f3.pack(fill="x", **pad)
        self.e_budget = self._entry(f3, "评估次数", "250", 6)
        self.e_time = self._entry(f3, "每次仿真(秒)", "8.0", 6)
        self.e_pop = self._entry(f3, "种群大小", "10", 5)
        self.e_seed = self._entry(f3, "随机种子", "1", 5)
        ttk.Label(f3, text="  可视化:").pack(side="left")
        self.cb_view = ttk.Combobox(f3, width=16, state="readonly",
                                    values=["每组都看", "每5组看1组", "只看刷新纪录", "不看(最快)"])
        self.cb_view.current(2)
        self.cb_view.pack(side="left", padx=4)
        self.v_rt = tk.BooleanVar(value=False)
        ttk.Checkbutton(f3, text="实时速度", variable=self.v_rt).pack(side="left", padx=6)

        # === 优化目标 & 搜索范围 ===
        f35 = ttk.LabelFrame(self.root, text=" 4. 优化目标 —— 决定“什么叫游得好” ")
        f35.pack(fill="x", **pad)
        r1 = ttk.Frame(f35); r1.pack(fill="x", pady=2)
        ttk.Label(r1, text="fitness =  速度  −").pack(side="left", padx=(8, 2))
        self.e_wyaw = ttk.Entry(r1, width=6); self.e_wyaw.insert(0, "0.3"); self.e_wyaw.pack(side="left")
        ttk.Label(r1, text="× 偏航率  −").pack(side="left", padx=2)
        self.e_wene = ttk.Entry(r1, width=6); self.e_wene.insert(0, "0.0"); self.e_wene.pack(side="left")
        ttk.Label(r1, text="× 平均功率(W)").pack(side="left", padx=2)
        ttk.Label(r1, text="     快捷:").pack(side="left", padx=(16, 4))
        for txt, wy, we in [("只求最快", "0", "0"), ("又快又直", "0.3", "0"),
                            ("严格直行", "1.0", "0"), ("省电优先", "0.3", "0.0003")]:
            ttk.Button(r1, text=txt, width=9,
                       command=lambda a=wy, b=we: self._set_w(a, b)).pack(side="left", padx=2)
        r2 = ttk.Frame(f35); r2.pack(fill="x", pady=2)
        ttk.Label(r2, text="搜索范围:   频率 f").pack(side="left", padx=(8, 2))
        self.e_f0 = ttk.Entry(r2, width=5); self.e_f0.insert(0, "0.2"); self.e_f0.pack(side="left")
        ttk.Label(r2, text="~").pack(side="left")
        self.e_f1 = ttk.Entry(r2, width=5); self.e_f1.insert(0, "2.5"); self.e_f1.pack(side="left")
        ttk.Label(r2, text="Hz     幅度 A 上限").pack(side="left", padx=(10, 2))
        self.e_amax = ttk.Entry(r2, width=5); self.e_amax.insert(0, "0.9"); self.e_amax.pack(side="left")
        ttk.Label(r2, text="rad  (调小可减少发散)").pack(side="left", padx=2)
        ttk.Label(r2, text="   偏置范围 ±").pack(side="left", padx=(10, 2))
        self.e_omax = ttk.Entry(r2, width=5); self.e_omax.insert(0, "0.5"); self.e_omax.pack(side="left")
        ttk.Label(r2, text="rad").pack(side="left")
        self.lb_fit = ttk.Label(f35, text="", foreground="#a0410d")
        self.lb_fit.pack(anchor="w", padx=10, pady=(2, 4))

        # === 按钮 ===
        f4 = ttk.Frame(self.root); f4.pack(fill="x", **pad)
        self.b_demo = ttk.Button(f4, text="▶ 看示例步态", command=self._demo)
        self.b_demo.pack(side="left", padx=4)
        self.b_run = ttk.Button(f4, text="🚀 开始寻优", command=self._start)
        self.b_run.pack(side="left", padx=4)
        self.b_stop = ttk.Button(f4, text="■ 停止", command=self._stop, state="disabled")
        self.b_stop.pack(side="left", padx=4)
        self.b_best = ttk.Button(f4, text="★ 看最优结果", command=self._play_best)
        self.b_best.pack(side="left", padx=4)
        self.lb_stat = ttk.Label(f4, text="就绪", foreground="#25507a", font=self.font_big)
        self.lb_stat.pack(side="left", padx=20)

        # === 中部: 左=参数表, 右=收敛曲线 ===
        mid = ttk.Frame(self.root); mid.pack(fill="both", expand=True, **pad)

        fl = ttk.LabelFrame(mid, text=" 当前最优步态参数 (每个电机的运动规律) ")
        fl.pack(side="left", fill="both", expand=True, padx=(0, 4))
        ttk.Label(fl, text="每个电机的角度 = 偏置 + 1次幅度·sin(2πf·t + 1次相位) + 2次幅度·sin(4πf·t + 2次相位)",
                  foreground="#555", wraplength=560).pack(anchor="w", padx=8, pady=(4, 0))
        ttk.Label(fl, text="偏置=摆动的中心角度(姿态) | 相位=起摆时刻(腿之间的相位差就是步态) | 2次谐波=脚蹼整流顺桨",
                  foreground="#888", wraplength=560).pack(anchor="w", padx=8, pady=(0, 4))
        cols = ("joint", "off", "a1", "p1", "a2", "p2")
        self.tv = ttk.Treeview(fl, columns=cols, show="headings", height=13)
        for c, t, w in zip(cols, ["电机", "偏置", "1次幅度", "1次相位", "2次幅度", "2次相位"],
                           [150, 80, 90, 90, 90, 90]):
            self.tv.heading(c, text=t); self.tv.column(c, width=w, anchor="center")
        self.tv.pack(fill="both", expand=True, padx=4, pady=4)
        self.lb_freq = ttk.Label(fl, text="频率: —", font=self.font_big, foreground="#0d4d2a")
        self.lb_freq.pack(anchor="w", padx=8, pady=(0, 6))

        fr = ttk.LabelFrame(mid, text=" 收敛情况 (曲线走平 + σ 变小 = 搜索充分) ")
        fr.pack(side="left", fill="both", expand=True, padx=(4, 0))
        self.cv = tk.Canvas(fr, bg="white", height=250)
        self.cv.pack(fill="both", expand=True, padx=6, pady=6)
        self.lb_conv = ttk.Label(fr, text="评估 0 | 最优 — | σ —", font=self.font_big,
                                 foreground="#0d4d2a")
        self.lb_conv.pack(anchor="w", padx=8, pady=(0, 6))

        # === 日志 ===
        f5 = ttk.LabelFrame(self.root, text=" 运行日志 ")
        f5.pack(fill="both", expand=True, **pad)
        self.txt = tk.Text(f5, height=9, font=self.font_mono)
        self.txt.pack(fill="both", expand=True, padx=4, pady=4)

    def _set_w(self, wy, we):
        self.e_wyaw.delete(0, "end"); self.e_wyaw.insert(0, wy)
        self.e_wene.delete(0, "end"); self.e_wene.insert(0, we)
        self._update_dim()

    def _fit_hint(self):
        try:
            wy, we = float(self.e_wyaw.get()), float(self.e_wene.get())
        except Exception:
            return "(权重填的不是数字)"
        if wy == 0 and we == 0:
            return "当前: 纯追求速度最快, 允许打转"
        if we > 0:
            return (f"当前: 兼顾省电 —— 偏航罚 {wy}, 功率罚 {we}/W "
                    f"(功率约 100~200 W 时罚 {we*150:.3f}, 与速度 0.05~0.1 同量级为宜)")
        if wy >= 1.0:
            return f"当前: 严格要求走直 (偏航罚 {wy}), 会牺牲一些速度"
        return f"当前: 速度与直行性折中 (偏航罚 {wy})"

    def _entry(self, parent, label, default, width):
        ttk.Label(parent, text=f"  {label}:").pack(side="left")
        e = ttk.Entry(parent, width=width); e.insert(0, default); e.pack(side="left", padx=2)
        return e

    # ---------------- 配置读写 ----------------
    def _read_cfg_text(self):
        return open(self.cfg_path, encoding="utf-8").read()

    def _load_cfg_into_ui(self):
        try:
            sys.path.insert(0, HERE)
            from simulate import load_cfg
            c = load_cfg(self.cfg_path)
            self.e_budget.delete(0, "end"); self.e_budget.insert(0, str(c.get("budget", 250)))
            self.e_time.delete(0, "end"); self.e_time.insert(0, str(c.get("sim_time", 8.0)))
            self.e_pop.delete(0, "end"); self.e_pop.insert(0, str(c.get("popsize", 10)))
            self.e_seed.delete(0, "end"); self.e_seed.insert(0, str(c.get("seed", 1)))
            self.e_wyaw.delete(0, "end"); self.e_wyaw.insert(0, str(c.get("w_yaw", 0.3)))
            self.e_wene.delete(0, "end"); self.e_wene.insert(0, str(c.get("w_energy", 0.0)))
            g = c.get("gait", {})
            fr = g.get("freq_range", [0.2, 2.5]); ar = g.get("amp_range", [0.0, 0.9])
            orr = g.get("offset_range", [-0.5, 0.5])
            for e, v in [(self.e_f0, fr[0]), (self.e_f1, fr[1]),
                         (self.e_amax, ar[1]), (self.e_omax, abs(orr[1]))]:
                e.delete(0, "end"); e.insert(0, str(v))
            m = c.get("model", "")
            if m in self.cb_model["values"]:
                self.cb_model.set(m)
        except Exception as e:
            self._log(f"[!] 读 config 失败: {e}")
        self._update_dim()

    def _make_run_cfg(self):
        """把界面上的设置写成临时配置 _ui_run.json (基于 config.json)."""
        from simulate import load_cfg, strip_jsonc
        c = load_cfg(self.cfg_path)
        c["budget"] = int(self.e_budget.get())
        c["sim_time"] = float(self.e_time.get())
        c["popsize"] = int(self.e_pop.get())
        c["seed"] = int(self.e_seed.get())
        c["w_yaw"] = float(self.e_wyaw.get())
        c["w_energy"] = float(self.e_wene.get())
        g = c.setdefault("gait", {})
        g["freq_range"] = [float(self.e_f0.get()), float(self.e_f1.get())]
        g["amp_range"] = [0.0, float(self.e_amax.get())]
        o = abs(float(self.e_omax.get()))
        g["offset_range"] = [-o, o]
        g["optimize"] = {"freq": self.v_freq.get(), "duty": self.v_duty.get(),
                         "amp": self.v_amp.get(),
                         "phase": self.v_phase.get(), "offset": self.v_off.get()}
        g["preset_phase"] = self.cb_preset.get().split()[0]
        json.dump(c, open("_ui_run.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        return "_ui_run.json"

    def _update_dim(self):
        try:
            from simulate import Swimmer, load_cfg
            p = self._make_run_cfg()
            sw = Swimmer(load_cfg(p)["model"], load_cfg(p))
            self.lb_dim.config(text=f"优化 {sw.gait.dim_opt} / {sw.gait.dim} 维")
            self.lb_model.config(text=f"执行器 {sw.model.nu} 个 | 刚体 {sw.model.nbody-1}")
            self.lb_fit.config(text=self._fit_hint())
        except Exception as e:
            self.lb_dim.config(text=f"(读取失败: {str(e)[:40]})")

    # ---------------- 模型 ----------------
    def _refresh_models(self):
        files = sorted(f.replace("\\", "/") for f in glob.glob("robots/*.xml"))
        self.cb_model["values"] = files
        if files and not self.cb_model.get():
            self.cb_model.current(0)

    def _switch_model(self):
        p = self.cb_model.get()
        r = subprocess.run([PY, "set_model.py", p], capture_output=True, text=True,
                           encoding="utf-8", errors="replace")
        self._log(r.stdout.strip() or r.stderr.strip())
        self._update_dim()

    def _import_urdf(self):
        f = filedialog.askopenfilename(title="选择 URDF", filetypes=[("URDF", "*.urdf"), ("所有文件", "*.*")])
        if not f: return
        name = os.path.splitext(os.path.basename(f))[0].lower()
        out = f"robots/{name}.xml"
        self._log(f"导入 {f} → {out} ...")
        r = subprocess.run([PY, "import_model.py", f, out, "--drive", "2.1,1.1",
                            "--armature", "3e-4", "--damping", "0.02", "--kp", "15", "--force", "3"],
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        self._log(r.stdout.strip() or r.stderr.strip())
        self._refresh_models()
        if os.path.exists(out):
            self.cb_model.set(out); self._switch_model()

    # ---------------- 运行 ----------------
    def _spawn(self, cmd, tag):
        if self.proc and self.proc.poll() is None:
            messagebox.showinfo("提示", "已有任务在跑, 先停止"); return
        self.txt.delete("1.0", "end"); self.hist.clear(); self._draw_curve()
        self._log("$ " + " ".join(cmd))
        self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, encoding="utf-8", errors="replace", bufsize=1)
        self.b_run.config(state="disabled"); self.b_stop.config(state="normal")
        self.lb_stat.config(text=tag)
        threading.Thread(target=self._reader, daemon=True).start()

    def _reader(self):
        for line in self.proc.stdout:
            self.q.put(line.rstrip())
        self.q.put(None)

    def _start(self):
        p = self._make_run_cfg()
        mode = self.cb_view.get()
        n = self.e_budget.get()
        if mode == "不看(最快)":
            cmd = [PY, "-u", "optimize.py", p, n]
        else:
            cmd = [PY, "-u", "optimize_view.py", p, n]
            if mode == "每5组看1组": cmd += ["--every", "5"]
            elif mode == "只看刷新纪录": cmd += ["--best"]
            if self.v_rt.get(): cmd += ["--rt"]
        self._spawn(cmd, "寻优中…")

    def _demo(self):
        self._spawn([PY, "-u", "view.py", self._make_run_cfg(), "--demo"], "播放示例步态")

    def _play_best(self):
        if not os.path.exists("results/best.json"):
            messagebox.showinfo("提示", "还没有 results/best.json, 先跑一次寻优"); return
        self._spawn([PY, "-u", "view.py", self._make_run_cfg(), "results/best.json"], "播放最优步态")

    def _stop(self):
        if self.proc and self.proc.poll() is None:
            self.proc.terminate()
        self.lb_stat.config(text="已停止")
        self.b_run.config(state="normal"); self.b_stop.config(state="disabled")

    # ---------------- 日志/曲线 ----------------
    def _log(self, s):
        if not s: return
        self.txt.insert("end", s + "\n"); self.txt.see("end")

    def _pump(self):
        try:
            while True:
                line = self.q.get_nowait()
                if line is None:
                    self.lb_stat.config(text="完成")
                    self.b_run.config(state="normal"); self.b_stop.config(state="disabled")
                    self._show_best()
                    continue
                self._log(line)
                m = re.search(r"#\s*(\d+).*?fitness\s*([-+0-9.]+)", line)
                if m:
                    ev, fit = int(m.group(1)), float(m.group(2))
                    best = max(fit, self.hist[-1][1] if self.hist else -9e9)
                    self.hist.append((ev, best)); self._draw_curve()
                m2 = re.search(r"evals\s+(\d+).*?最优\s*([-+0-9.]+)|已评估\s*(\d+).*?最优\s*([-+0-9.]+)", line)
                ms = re.search(r"sigma\s+([0-9.]+)", line)
                if ms or m2:
                    sig = f"{float(ms.group(1)):.4f}" if ms else "—"
                    ev = self.hist[-1][0] if self.hist else 0
                    bf = f"{self.hist[-1][1]:+.4f}" if self.hist else "—"
                    self.lb_conv.config(text=f"评估 {ev} | 最优 {bf} | σ {sig}")
                if "★" in line or "新纪录" in line:
                    self._show_best()
        except queue.Empty:
            pass
        self.root.after(250, self._pump)

    def _draw_curve(self):
        c = self.cv; c.delete("all")
        W = c.winfo_width() or 500; H = c.winfo_height() or 250
        if len(self.hist) < 2:
            c.create_text(W/2, H/2, text="等待数据…", fill="#888"); return
        pts_ok = [(x, y) for x, y in self.hist if y > -1.0]      # 丢掉发散(-1000)的点
        if len(pts_ok) < 2:
            c.create_text(W/2, H/2, text="等待有效数据…(前几次可能发散)", fill="#888")
            return
        xs = [p[0] for p in pts_ok]; ys = [p[1] for p in pts_ok]
        x0, x1 = min(xs), max(xs); y0, y1 = min(ys), max(ys)
        span = max(y1 - y0, 1e-3); y0 -= span*0.1; y1 += span*0.1
        mL, mR, mT, mB = 46, 12, 12, 26
        def px(x): return mL + (x-x0)/(x1-x0+1e-9)*(W-mL-mR)
        def py(y): return H-mB - (y-y0)/(y1-y0)*(H-mT-mB)
        c.create_line(mL, H-mB, W-mR, H-mB, fill="#999")
        c.create_line(mL, mT, mL, H-mB, fill="#999")
        pts = []
        for x, y in pts_ok: pts += [px(x), py(y)]
        c.create_line(*pts, fill="#c81e3c", width=3)
        c.create_oval(px(xs[-1])-4, py(ys[-1])-4, px(xs[-1])+4, py(ys[-1])+4,
                      fill="#c81e3c", outline="")
        c.create_text(W/2, H-8, text="评估次数 →", fill="#555", font=("", 9))
        c.create_text(mL-8, mT+8, text=f"{y1:+.3f}", anchor="e", fill="#555", font=("", 9))
        c.create_text(mL-8, H-mB, text=f"{y0:+.3f}", anchor="e", fill="#555", font=("", 9))
        c.create_text(W-mR, H-mB+14, text=f"{x1}", anchor="e", fill="#555", font=("", 9))
        c.create_text(mL+6, mT+8, text="fitness (越高越好)", anchor="w", fill="#999", font=("", 9))

    def _show_best(self):
        try:
            from simulate import Swimmer, load_cfg
            b = json.load(open("results/best.json", encoding="utf-8"))
            c = load_cfg("_ui_run.json" if os.path.exists("_ui_run.json") else self.cfg_path)
            sw = Swimmer(c["model"], c)
            p = sw.gait.decode(sw.gait.expand(b["x"]))
            self.lb_freq.config(text=f"频率: {p['freq']:.3f} Hz    "
                                     f"速度: {b.get('speed',0):+.4f} m/s    偏航: {b.get('yaw',0):.1f}°")
            for it in self.tv.get_children(): self.tv.delete(it)
            import numpy as np
            for i, nm in enumerate(sw.gait.names):
                off = p["offs"][sw.gait._off_slot[i]]
                a1 = p["A"][0][sw.gait._amp_slot[i]]; p1 = np.degrees(p["PH"][0][i])
                a2 = p["A"][1][sw.gait._amp_slot[i]] if sw.gait.H > 1 else 0
                p2 = np.degrees(p["PH"][1][i]) if sw.gait.H > 1 else 0
                self.tv.insert("", "end", values=(nm, f"{off:+.3f}", f"{a1:.3f}", f"{p1:6.1f}°",
                                                  f"{a2:.3f}", f"{p2:6.1f}°"))
        except Exception as e:
            self._log(f"[!] 读最优失败: {e}")


if __name__ == "__main__":
    root = tk.Tk()
    try: root.call("tk", "scaling", 1.15)
    except Exception: pass
    App(root)
    root.mainloop()
