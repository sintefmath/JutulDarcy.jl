<template>
  <span
    class="badge"
    :style="{ backgroundColor: color }"
  >
    <a :href="`../../examples/overview/example_overview#${text}`">{{ text }}</a> 
  </span>
  &nbsp;
</template>

<script setup>
defineProps({
  color: {
    type: String,
    required: true
  },
  text: {
    type: String,
    required: true
  }
});
</script>

<style scoped>
.badge {
  display: inline-flex;
  align-items: center;
  gap: 0.25em;

  padding: 0.2em 0.5em;
  border-radius: 999px;

  font-size: 0.85em;
  font-weight: 500;
  line-height: 1;

  user-select: none;
  white-space: nowrap;

  text-decoration: none;   /* remove link underline */
  color: inherit;          /* prevent default link color */

  background-color: var(--badge-bg, rgba(0, 0, 0, 0.08));
  transition: background-color 0.15s ease, transform 0.05s ease;
}

/* hover / focus states */
.badge:hover {
  background-color: var(--badge-bg-hover, rgba(0, 0, 0, 0.12));
  text-decoration: none;
}

.badge:active {
  transform: scale(1.0);
}

/* optional: if it's an <a> element */
a.badge {
  color: inherit;
  text-decoration: none;
}

.badge,
.badge * {
  text-decoration: none !important;
  color: inherit !important;
}
</style>

<!-- Example usage -->
<!--
<ExampleTag color="#3498db" text="Info" />
<ExampleTag color="green" text="Success" />
-->
