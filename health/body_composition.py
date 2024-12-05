import pandas as pd
import seaborn as sns
import datetime as dt
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 8})
sns.set_theme(style="darkgrid")
sns.set_context("notebook", font_scale=0.8)
colors = ["#a84432", "#ebe698"]

# Read the Excel file
file_path = "./body_composition.xlsx"
df = pd.read_excel(file_path)
df = df.transpose()
df.columns = df.iloc[0, :]
df = df.iloc[1:, :]
df.index = [dt.datetime.strptime(t, "%b %Y") for t in df.index.values]
df.pop("Scan type")


fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(18, 8), sharex=True, sharey=True)

# Flatten the axes array for easy iteration
axes = axes.flatten()

# Plot each body part in a subplot
for i, body_part in enumerate(
    [
        "Left arm",
        "Right arm",
        "Toraxic spine",
        "Left ribs",
        "Right ribs",
        "Lumbar spine",
        "Left leg",
        "Right leg",
        "Pelvis",
    ]
):
    df[[f"{body_part} muscle (kg)", f"{body_part} fat (kg)"]].plot(
        kind="area",
        stacked=True,
        alpha=0.5,
        color=colors,
        ax=axes[i],
        legend=False,
    )
    axes[i].set_title(body_part)

# Set common labels
fig.text(0.001, 0.5, "Mass (kg)", va="center", rotation="vertical")
plt.tight_layout()


########################################################################
########################################################################

fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(14, 6))

# # # Flatten the axes array for easy iteration
axes = axes.flatten()

# Plot each body part in a subplot
parameters = ["Visceral fat area (cm2)", "Subcutaneous fat area (cm2)"]
df[parameters].plot(
    kind="area",
    stacked=True,
    alpha=0.5,
    ax=axes[0],
    color=[
        "#3f52b0",
        "#98c4eb",
    ],
)

df["Weight (kg)"].plot(
    kind="area", alpha=0.5, ax=axes[1], color=["#3f52b0"], legend="lower left"
)

df["% Fat"].plot(
    kind="area", alpha=0.5, ax=axes[2], color=["#f7d67c"], legend="lower left"
)
df["FMI (kg/m2)"].plot(
    kind="area", alpha=0.5, ax=axes[3], color=["#f7d67c"], legend="lower left"
)
df["Waist/Hip ratio"].plot(
    kind="area", alpha=0.5, ax=axes[4], color=["#4da39c"], legend="lower left"
)
axes[4].plot([df.index.values[0], df.index.values[-1]], [1, 1], "k--")
df["Lean mass (kg)"].plot(
    kind="area", alpha=0.5, ax=axes[5], color=["#a84432"], legend="lower left"
)

plt.tight_layout()

plt.show()
