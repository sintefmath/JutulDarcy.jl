<template>
  <span
    class="badge"
    :style="{ backgroundColor: color }"
  >
    <a href="man/basics/input_files">{{ text }}</a>
  </span>
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
  display: inline-block;
  padding: 0.3em 0.8em;
  border-radius: 999px;
  color: #fff;
  font-weight: 500;
  font-size: 1em;
  user-select: none;
}
</style>

<!-- Example usage -->
<!--
<ExampleTag color="#3498db" text="Info" />
<ExampleTag color="green" text="Success" />
-->
